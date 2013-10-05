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
#include <boost/algorithm/string/predicate.hpp>

struct VertexInfo
{
  std::string name;
};

typedef boost::property<boost::edge_weight_t, int> EdgeWeight; 

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexInfo, EdgeWeight> Graph;
typedef Graph::vertex_descriptor Vertex; 

// Структура описывает фильтр для вершины.
// Под фильтр попадает вершина, если ее имя начинается с Prefix
// и число соединений равно Connections.
struct VertexFilter
{
  VertexFilter(const std::string& prefix, int connections)
  :Prefix(prefix),Connections(connections){}
  
  VertexFilter(const std::string& filter):Connections(0)
  {
    size_t numPos = filter.find_first_of("1234567890");
    
    Prefix = filter.substr(0, numPos);
    
    if (numPos != std::string::npos)
    {
      Connections = atoi(filter.substr(numPos).c_str());
    }
  }
  
  bool match(const Vertex& v, const Graph& graph)
  {
    int vConnections = boost::in_degree(v,  graph); // число ребер для вершины v
    
    return ((Connections == 0 || vConnections == Connections) 
    && (Prefix.empty() || boost::starts_with(graph[v].name, Prefix)));
  }
  
  bool operator==(const VertexFilter& other) const
  {
    return Prefix == other.Prefix && Connections == other.Connections;
  }
  
  friend std::ostream& operator<<(std::ostream& stream, const VertexFilter& filter)
  {
    stream << filter.Prefix;
    if (filter.Connections>0)
    {
      stream << filter.Connections;
    }
    return stream;
  }
  
  std::string Prefix;
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
  void WriteGraphviz(std::ostream& out) const
  {
    boost::write_graphviz(out, _graph, boost::make_label_writer(boost::get(&VertexInfo::name, _graph)));
  }
  
  size_t CountVertices(const VertexFilter& filter) const
  {
    std::vector<Vertex> vertices = findAll(filter);
    return vertices.size();
  }
  
  // Cчитает число связей между узлами, подходящими под фильтр.
  // Подходит для расчета Cn-Cn. Не походит для Cn-Cm, где n != m
  size_t CountConnections(const VertexFilter& filter) const
  {
    std::vector<Vertex> vertices = findAll(filter);
        
    size_t result = 0;
    if (!vertices.empty())
    {
      for (size_t i=0; i<vertices.size()-1; ++i)
      {
	for (size_t j=i+1; j<vertices.size(); ++j)
	{
	  result += Distance(vertices[i], vertices[j]);
	}
      }
    }
    
    return result;
  }
  
  // Cчитает число связей от атомами, подходящими под фильтр from.
  // до атомов, подходящих под фильтр to
  // Подходит как для расчета Cn-Cm, где n != m, так и для Cn-Cn, Cn-OH и т.п.
  size_t CountConnections(const VertexFilter& from, const VertexFilter& to) const
  {
    if (from == to)
    {
      return CountConnections(from);
    }
        
    std::vector<Vertex> fromVertices = findAll(from);
    std::vector<Vertex> toVertices = findAll(to);
    
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
  
  // Ищет все узлы, назвиния которых начинаются с prefix
  // и имеют число связей равное connections
  std::vector<Vertex> findAll(const VertexFilter& filter) const
  {
    std::vector<Vertex> result;
    
    Graph::vertex_iterator vBegin, vEnd;
    boost::tie(vBegin, vEnd) = boost::vertices(_graph);
   
    std::copy_if(vBegin, vEnd, std::back_inserter(result), 
		 std::bind(&VertexFilter::match, filter, std::placeholders::_1, _graph));
    
    return result;
  }
  
  // Считает число ребер между двумя вершинами
  // Магия boost graph library
  size_t Distance(Vertex v1, Vertex v2) const
  {
    std::vector<Vertex> predecessors;
    predecessors.resize(num_vertices(_graph));
    
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

// Проводит группу вычислений с указанной молекулой
// Вычисления счтитываются из файла.
// Формат файла -- одно вычисление на одной строке.
// Формат описания вычисления
//  -- либо один фильтр вида <Prefix><Connections>, например C1,C4 или OH
//  -- либо два фильтра вида <Prefix><Connections>-<Prefix><Connections>,
//	например С1-С1, С1-С3, С2-OH
// Для одиночных фильтров будет посчитано общее число таких вершин.
// Для двойных -- число связей между вершинами, подходящими под первый и второй фильтр
class Calculator
{
public:
  Calculator(std::ifstream& file)
  {
    while(file.peek() != EOF)
    {
      std::string calculation;
      std::getline(file, calculation);
      
      if (calculation.empty())
	continue;
      
      size_t splitPos = calculation.find("-");
      if (splitPos != std::string::npos)
      {
	VertexFilter from(calculation.substr(0, splitPos));
	VertexFilter to(calculation.substr(splitPos+1));
	
	_connectionCalcs.push_back(std::pair<VertexFilter, VertexFilter>(from,to));
      }
      else
      {
	_vertexCalcs.push_back(VertexFilter(calculation));
      }
    }
  }
  
  // Печатает заголовок таблицы с перечислеием всех вычислений.
  // Сначала идут вычисления для числа вершин,
  // потом -- для числа связей.
  // Порядок внутри каждой группы такой же, как во входном файле.
  void printHeaders(std::ostream& stream, char delimiter = '\t') const
  {
    for (auto it = _vertexCalcs.cbegin(); it != _vertexCalcs.cend(); ++it)
    {
      stream << delimiter << *it;
    }
    
    for (auto it = _connectionCalcs.cbegin(); it != _connectionCalcs.cend(); ++it)
    {
      stream << delimiter << it->first;
      stream << "-" << it->second;
    }
  }
  
  // Возвращает набор результатов вычислений.
  // Сначала идут вычисления для числа вершин,
  // потом -- для числа связей.
  // Порядок внутри каждой группы такой же, как во входном файле.
  std::vector<size_t> calculate(const Molecule& m) const
  {
    std::vector<size_t> result;
    
    for (auto it = _vertexCalcs.cbegin(); it != _vertexCalcs.cend(); ++it)
    {
      result.push_back(m.CountVertices(*it));
    }
    
    for (auto it = _connectionCalcs.cbegin(); it != _connectionCalcs.cend(); ++it)
    {
      result.push_back(m.CountConnections(it->first, it->second));
    }
    
    return result;
  }
  
private:
  std::list<VertexFilter> _vertexCalcs;
  std::list<std::pair<VertexFilter, VertexFilter>> _connectionCalcs;
};

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::cout<<"Usage: prog formula_file calc_file "<<std::endl;
    std::cout<<"formula_file -- file with formulas, one formula on each line"<<std::endl;
    std::cout<<"calc_file -- file with parameters to calculate, one parameter on line"<<std::endl;
    std::cout<<"Program creates file filename.csv -- result of calculations"<<std::endl;
    std::cout<<"Program creates files filename_i.dot and filename_i.png -- images of graphs for each formula"<<std::endl;
    return 0;
  }
  
  std::string filePrefix = argv[1];
  
  std::ifstream inFile(argv[1]);
  
  std::ifstream calcFile(argv[2]);
  
  Calculator calculator(calcFile);
  
  std::ofstream outFile(filePrefix + "_out.csv");
    
  // Заголовок таблицы в CSV
  outFile << "N\tFormula";
  calculator.printHeaders(outFile);
  outFile<<std::endl;
  
  int formulaCount = 0;
  
  while (inFile.peek() != EOF)
  {
    std::string formula;
    std::getline(inFile, formula);
    
    if (formula.empty()) continue;
    
    formulaCount++;
    
    Molecule m(formula);
    
    std::string dotFileName = (boost::format("%1%_%2%.dot") % filePrefix % formulaCount).str();
    std::string pngFileName = (boost::format("%1%_%2%.png") % filePrefix % formulaCount).str();
    
    std::ofstream dotFile(dotFileName);
    m.WriteGraphviz(dotFile);
    dotFile.close();
    
    std::system((boost::format("dot -Tpng -o %1% %2%") % pngFileName % dotFileName).str().c_str());
    
    outFile << formulaCount << '\t' << formula;
    std::vector<size_t> results = calculator.calculate(m);
    
    for (auto it = results.cbegin(); it != results.cend(); ++it)
    {
      outFile << '\t' << *it;
    }
    outFile << std::endl;
  }
  
  inFile.close();
  calcFile.close();
  outFile.close();
  
  return 0;
}

