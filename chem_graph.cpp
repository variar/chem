/*
 Copyright (c) 2013, Anton Filimonov <anton.filimonov@gmail.com>
 A ll rights reserved.                                *
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met: 
 
 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer. 
 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution. 
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 Needs boost with graph library and graphviz installed
 compile: g++ -std=c++11 chem_graph.cpp
*/

#include <iostream>
#include <fstream>
#include <stack>
#include <list>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/optional.hpp>
#include <boost/format.hpp>

struct VertexInfo
{
  std::string name;
};

typedef boost::property<boost::edge_weight_t, int> EdgeWeight; 

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexInfo, EdgeWeight> Graph;
typedef Graph::vertex_descriptor Vertex; 

struct VertexFilter
{
  VertexFilter(const std::string& elements, int connections = 1)
  :Elements(elements),Connections(connections){}
  
  std::string Elements;
  int Connections;
};

class Molecule
{
public:
  Molecule(const std::string& formula, const std::string& baseElements = "CON" )
  :_formula(formula),_baseElements(baseElements)
  {
    BuildGraph();
  }
  
  // Записывает в выходной поток граф в формате graphviz (на языке dot)
  // Потом можно построить картинку с помощью программы dot
  void WriteGraphviz(std::ostream& out)
  {
    boost::write_graphviz(out, _graph, boost::make_label_writer(boost::get(&VertexInfo::name, _graph)));
  }
  
  size_t VertexCount()
  {
    return num_vertices(_graph);
  }
  
  // Считает все связи между атомами, перечисленными в строке elements,
  // у которых число связей равно connections.
  // Для CH3 или OH connections = 1, для CH2 -- 2 и т.д.
  size_t CountAllConnections(int connections, const std::string& elements)
  {
    std::vector<Vertex> vertices = findAll(connections, elements);
    
    size_t result = 0;
    for (size_t i=0; i<vertices.size()-1; ++i)
    {
      for (size_t j=i+1; j<vertices.size(); ++j)
      {
	result += Distance(vertices[i], vertices[j]);
      }
    }
    
    return result;
  }
  
  size_t CountConnections(const VertexFilter& from, const VertexFilter& to)
  {
    std::vector<Vertex> fromVertices = findAll(from.Connections, from.Elements);
    std::vector<Vertex> toVertices = findAll(to.Connections, to.Elements);
    
    size_t result = 0;
    for (size_t i=0; i<fromVertices.size(); ++i)
    {
      for (size_t j=0; j<toVertices.size(); ++j)
      {
	result += Distance(fromVertices[i], toVertices[j]);
      }
    }
    
    return result;
  }
  
private:
  
  // Строит внутреннее представление графа, разбирая строку формулы.
  // Символам C и O в строке будут соответствовать узлы графа.
  // Если в формуле встречается гурппа в (), то она присоединяется к 
  // узлу, стоящему перед "(". Следующий за ")"  узел будет присоединен
  // к узлу, стоящему перед "(".
  void BuildGraph()
  {
    boost::optional<Vertex> lastVertex; // последняя вершина из цепочки, в начале тут пусто.
    std::stack<Vertex> subgroupVertex; // стек вершин, от которых ответвляются группы (нужен для обработки вложенных () ) 
    
    // Cюда собирается название вершины. Берется часть строки от элемента C или O
    // и до следующего элемента С или О (или до символов () )
    std::string vertexName; 
    
    // вспомогательная функция, которая устанавливает назание для вершины,
    // если для нее еще не задано название. Нужна для корректной обработки подгрупп в ().
    auto setVertexName = [] (Graph& g, const Vertex& v, std::string name ) -> bool
    {
      bool nameSet = false;
      if (g[v].name.empty())
      {
	g[v].name = name;
	nameSet = true;
      }
      return nameSet;
    };
    
    // Цикл по каждому символу в строке формулы
    for (int i=0; i<_formula.length(); ++i)
    {
      // Для базовых атомов (C и O) создаем новые узлы
      if (_baseElements.find(_formula[i]) != std::string::npos)
      {
	//found new vertex
	Vertex newVertex = boost::add_vertex( _graph);
	
	// если у нас есть предыдущий узел, то добавляем ребро к новому
	if (lastVertex) 
	{
	  boost::add_edge(*lastVertex, newVertex, 1, _graph);
	  
	  // Т.к. мы нашли новый узел, то в vertexName собралось название предыдущей вершины
	  // Если у нее еще не установлено название, то устанавливаем и очищаем строку
	  // для сбора следующего названия.
	  if (setVertexName(_graph, *lastVertex, vertexName))
	  {
	    vertexName.clear();
	  }
	}
	
	// теперь новая вершина в хвосте цепочки
	lastVertex.reset(newVertex);
      }
      else if (_formula[i] == '(')
      {
	// начинается новая подгруппа -- ответвление
	// можно установить название для предыдущего узла
	if (setVertexName(_graph, *lastVertex, vertexName))
	{
	  vertexName.clear();
	}
	
	// мы запоминаем узел, от которого началось ответвление,
	// чтобы продолжить цепочку с него, когда ответвление закончится
	subgroupVertex.push(*lastVertex);
      }
      else if (_formula[i] == ')')
      {
	// группа кончилась, можно установить название ее последнего узла
	if (setVertexName(_graph, *lastVertex, vertexName))
	{
	  vertexName.clear();
	}
	
	// мы достаем из стека узел, от которого началось ветвление,
	// чтобы продолжить строить цепочку с него.
	lastVertex.reset(subgroupVertex.top());
	subgroupVertex.pop();
      }
      
      if (_formula[i] != '(' && _formula[i] != ')')
	vertexName += _formula[i];
    }
    
    _graph[*lastVertex].name = vertexName;
  }
  
  // Ищет все узлы, которые соответствуют элементам из elements
  // и имеют число связей равное connections
  std::vector<Vertex> findAll(int connections, const std::string& elements)
  {
    // вспомогательная функция - предикат, которая позволяет отбирать узлы из графа
    auto predicate = [this, connections, &elements](const Vertex& v) -> bool
    {
      int vConnections = boost::in_degree(v,  _graph); // число ребер для вершины v
      char vElement = _graph[v].name.empty() ? '\0' : _graph[v].name[0];
      
      return (vConnections == connections && (elements.empty() || elements.find(vElement) != std::string::npos));
    };
    
    std::vector<Vertex> result;
    
    Graph::vertex_iterator vBegin, vEnd;
    boost::tie(vBegin, vEnd) = boost::vertices(_graph);
   
    std::copy_if(vBegin, vEnd, std::back_inserter(result), predicate);
    
    return result;
  }
  
  // Считает число ребер между двумя вершинами
  // Магия boost graph library
  size_t Distance(Vertex v1, Vertex v2)
  {
    std::vector<Vertex> predecessors;
    predecessors.resize(VertexCount());
    
    boost::breadth_first_search(_graph,
				v1,
				boost::visitor(
				  boost::make_bfs_visitor(
				    boost::record_predecessors(&predecessors[0],boost::on_tree_edge())
				  )
				));
    
    size_t d = 0;
    for (Vertex v = v2; v != v1; v = predecessors[v])
    {
      d++;
    }
    
    return d;
  }
  
  std::string _formula;
  std::string _baseElements;
  Graph _graph;
};

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cout<<"Usage: prog filename"<<std::endl;
    std::cout<<"filename -- file with formulas, one formula on each line"<<std::endl;
    std::cout<<"Program creates file filename.csv -- result of calculations"<<std::endl;
    std::cout<<"Program creates files filename_i.dot and filename_i.png -- images of graphs for each formula"<<std::endl;
    return 0;
  }
  
  std::string filePrefix = argv[1];
  
  std::ifstream inFile(argv[1]);
  
  std::ofstream outFile(filePrefix + "_out.csv");
  
  
  // Заголовок таблицы в CSV
  outFile << "N\tFormula\tC1-OH\tC2-OH\tC3-OH\tC4-OH"<<std::endl;
  
  int formulaCount = 0;
  
  while (inFile.peek() != EOF)
  {
    std::string formula;
    std::getline(inFile, formula);
    formulaCount++;
    
    Molecule m(formula);

    outFile << formulaCount << "\t" << formula;
    
    VertexFilter filterOH("O",1);
    for (size_t i = 1; i<5; ++i)
    {
      outFile << "\t" << m.CountConnections(VertexFilter("C",i), filterOH);
    }
    outFile << std::endl;
    
    std::string dotFileName = (boost::format("%1%_%2%.dot") % filePrefix % formulaCount).str();
    std::string pngFileName = (boost::format("%1%_%2%.png") % filePrefix % formulaCount).str();
    
    std::ofstream dotFile(dotFileName);
    m.WriteGraphviz(dotFile);
    dotFile.close();
    
    std::system((boost::format("dot -Tpng -o %1% %2%") % pngFileName % dotFileName).str().c_str());
  }
  
  inFile.close();
  outFile.close();
  
  return 0;
}

