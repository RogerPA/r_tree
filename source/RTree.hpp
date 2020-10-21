// Copyright

#ifndef SOURCE_RTREE_HPP_
#define SOURCE_RTREE_HPP_

#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

#include "Rectangle.hpp"
using namespace std;

template <size_t N, typename ElemType, size_t M, size_t m = M / 2>
class RTree {
 public:
  struct Node;

  struct SpatialObject {
    Rectangle<N> box;
    ElemType identifier;
    std::shared_ptr<Node> child_pointer;
  };

  struct Node {
    typedef SpatialObject *iterator;
    typedef const SpatialObject *const_iterator;

    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;

    SpatialObject &operator[](size_t index);
    SpatialObject operator[](size_t index) const;

    bool is_leaf();

    std::shared_ptr<Node> insert(const SpatialObject &new_entry);
	
    SpatialObject entry[M];
    size_t size = 0;
  };

  RTree();
  virtual ~RTree();
  size_t dimension() const;
  size_t size() const;
  bool empty() const;

  void insert(const Rectangle<N> &box, const ElemType &value);
  std::shared_ptr<Node> choose_leaf(const std::shared_ptr<Node> &current_node,
                                    const Rectangle<N> &box,
                                    const ElemType &value);

  std::shared_ptr<Node> choose_node(const std::shared_ptr<Node> &current_node,
                                    const Rectangle<N> &box,
                                    SpatialObject *&entry);

  std::shared_ptr<Node> adjust_tree(const std::shared_ptr<Node> &parent,
                                    const std::shared_ptr<Node> &left,
                                    const std::shared_ptr<Node> &right,
                                    SpatialObject *entry);

  // TODO(ADE): Implement the details of all this functions
  std::vector<ElemType> &operator[](const Rectangle<N> &box);
  std::vector<ElemType> &at(const Rectangle<N> &box);
  const std::vector<ElemType> &at(const Rectangle<N> &box) const;
  void buscar(Rectangle<N> &rectangulo,shared_ptr<Node> nodo_actual=nullptr);
  // std::vector<ElemType> kNNValue(const Rectangle<N> &box, size_t k) const;

  // private:
  std::shared_ptr<Node> root_pointer_;
  size_t contador;
  vector<ElemType> historial;
};

/** Node R-tree struct implementation details*/
template <size_t N, typename ElemType, size_t M, size_t m>
typename RTree<N, ElemType, M, m>::Node::iterator
RTree<N, ElemType, M, m>::Node::begin() {
  return entry;
}

template <size_t N, typename ElemType, size_t M, size_t m>
typename RTree<N, ElemType, M, m>::Node::iterator
RTree<N, ElemType, M, m>::Node::end() {
  return entry + size;
}

template <size_t N, typename ElemType, size_t M, size_t m>
typename RTree<N, ElemType, M, m>::Node::const_iterator
RTree<N, ElemType, M, m>::Node::begin() const {
  return entry;
}

template <size_t N, typename ElemType, size_t M, size_t m>
typename RTree<N, ElemType, M, m>::Node::const_iterator
RTree<N, ElemType, M, m>::Node::end() const {
  return entry + size;
}

template <size_t N, typename ElemType, size_t M, size_t m>
typename RTree<N, ElemType, M, m>::SpatialObject
    &RTree<N, ElemType, M, m>::Node::operator[](size_t index) {
  return entry[index];
}

template <size_t N, typename ElemType, size_t M, size_t m>
typename RTree<N, ElemType, M, m>::SpatialObject
    RTree<N, ElemType, M, m>::Node::operator[](size_t index) const {
  return entry[index];
}

template <size_t N, typename ElemType, size_t M, size_t m>
bool RTree<N, ElemType, M, m>::Node::is_leaf() {
  if (size && entry[0].child_pointer) {
    return false;
  }
  return true;
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::Node::insert(const SpatialObject &new_entry) {
  if (size < M) {
    entry[size++] = new_entry;
    return nullptr;
  }
  /*
  ////////////////////////////////////////////cuadratico
  dividir un conjunto de M+1 entradas ded indice en 2 entradas
  //usa pickseeds para elegir 2 entradas que ser�n los primeros elementos de los entradas. 
  //asignar cada una a un grupo
  
  //Si todas las entradas han sido asignadas
    //para
  //si un grupo tiene pocas entradas que todas las demas asignarlo en el orden m  para que tenga el numero minimo 
    //asignarlo y para
  //usa Pick_Next para elegir la siguiente entrada a asignar.
  //a��dalo al grupo cuyo rectangulo de cobertura tendra que ser ampliado lo menos posible
  para acomodarlo  
  //Resuelva los v�nculos a�adiendo la entrada al grupo con un area mas pequena, 
  luego al de menos entradas, y luego a cualquiera de los dos 
  */
  vector<SpatialObject> entradas(M+1);
  entradas[M]=new_entry;
  size=0;
  int entrada1,entrada2,n_id;
  
  double area;
  double incre_area_1,incre_area_2;
  
  Rectangle<N> rectangulo,rec1,rec2;
  for(int i=0;i<M;++i){
	for(int j=i+1;j<M;++j){  
		rectangulo = entradas[i].box;
		rectangulo.adjust(entradas[j].box);
		area = rectangulo.get_area() - entradas[i].box.get_area() - entradas[j].box.get_area();
		if (area > 0) {
			entrada1 = i;
			entrada2 = j;
		}
	}
  }
  insert(entradas[entrada1]);
  
  shared_ptr<Node> nodo_nuevo;
  nodo_nuevo->insert(entradas[entrada2]);
  
  Rectangle<N> grupo1=entradas[entrada1].box;
  Rectangle<N> grupo2=entradas[entrada2].box;
  
  entradas.erase(entradas.begin()+entrada1,entradas.begin()+entrada2-1);
  
  while(entradas.size()!=0){
	if(size+entradas.size()==m){
		for(int i=0;i<entradas.size();++i){
			insert(entradas[i]);
		}
		entradas.clear();
	}
	else if (nodo_nuevo->size + entradas.size() ==m){
		for(int i=0;i<entradas.size();++i){
			nodo_nuevo->insert(entradas[i]);
		}
		entradas.clear();
	}
	else{
		double diferencia=0;
		Rectangle<N> g1,g2;
		double area_g1,area_g2,diferencia_area;
		for (int i = 0; i < entradas.size(); ++i) {
			g1 = grupo1;
			g1.adjust(entradas.at(i).box);
			area_g1 = g1.get_area() - grupo1.get_area();
			
			g2 = grupo2;
			g2.adjust(entradas[i].box);
			area_g2 = g2.get_area() - grupo2.get_area();
			
			diferencia_area = abs(area_g1 - area_g2);
			
			if (diferencia_area >= diferencia) {
				diferencia = diferencia_area;
				n_id = i;
			}
		}
		rec1=grupo1;
		rec1.adjust(entradas[n_id].box);
		incre_area_1=rec1.get_area() - grupo1.get_area();
		
		rec2=grupo2;
		rec2.adjust(entradas[n_id].box);
		incre_area_2=rec2.get_area() - grupo2.get_area();
		if(incre_area_1 < incre_area_2 || grupo1.get_area()<grupo2.get_area()){
			grupo1.adjust(entradas[n_id].box);
			insert(entradas[n_id]);
		}
		else{
			grupo2.adjust(entradas[n_id].box);
			nodo_nuevo->insert(entradas[n_id]);
		}
		entradas.erase(entradas.begin()+n_id);
	  }
    }
	return nodo_nuevo;
} 

/** R-Tree class implementation details */
///
template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::RTree() : root_pointer_(new Node) {
	contador=0;
}

// TODO(ADE):
template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::~RTree() {}
///
template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::dimension() const {
  return N;
}
///
template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::size() const {
  return contador;
}
///
template <size_t N, typename ElemType, size_t M, size_t m>
bool RTree<N, ElemType, M, m>::empty() const {
  return size()==0;
}
///
template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::insert(const Rectangle<N> &box,
                                      const ElemType &value) {
  contador+=1;
  std::shared_ptr<Node> splitted_node = choose_leaf(root_pointer_, box, value);
  if (!splitted_node) {
    return;
  }
  shared_ptr<Node> root_temporal=make_shared<Node>();
  root_temporal->entry[0].child_pointer = root_pointer_;
  root_temporal->size+=1;
  adjust_tree(root_temporal, root_pointer_, splitted_node, &root_temporal->entry[0]);
  
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::choose_leaf(const std::shared_ptr<Node> &current_node,
                                      const Rectangle<N> &box,
                                      const ElemType &value) {
  if (!current_node->is_leaf()) {
    SpatialObject *entry;
    std::shared_ptr<Node> next_node = choose_node(current_node, box, entry);
    std::shared_ptr<Node> splitted_node = choose_leaf(next_node, box, value);
    return adjust_tree(current_node, next_node, splitted_node, entry);
  }
  SpatialObject new_entry;
  new_entry.box = box;
  new_entry.identifier = value;
  return current_node->insert(new_entry);
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::choose_node(const std::shared_ptr<Node> &current_node,
                                      const Rectangle<N> &box,
                                      SpatialObject *&entry) {
  float minimum_area = (*current_node)[0].box.get_area();

  Rectangle<N> enlarged_box = (*current_node)[0].box;
  enlarged_box.adjust(box);
  float minimum_enlargement = enlarged_box.get_area() - minimum_area;

  float enlargement, area;
  std::shared_ptr<Node> node = (*current_node)[0].child_pointer;

  entry = &(*current_node)[0];
  for (SpatialObject &current_entry : *current_node) {
    area = current_entry.box.get_area();

    enlarged_box = current_entry.box;
    enlarged_box.adjust(box);
    enlargement = enlarged_box.get_area() - area;

    if (enlargement < minimum_enlargement ||
        (enlargement == minimum_enlargement && area < minimum_area)) {
      minimum_enlargement = enlargement;
      minimum_area = area;
      node = current_entry.child_pointer;
      entry = &current_entry;
    }
  }

  return node;
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::adjust_tree(const std::shared_ptr<Node> &parent,
                                      const std::shared_ptr<Node> &left,
                                      const std::shared_ptr<Node> &right,
                                      SpatialObject *entry) {
  entry->box.reset();
  for (SpatialObject current_entry : *left) {
    entry->box.adjust(current_entry.box);
  }
  if (!right) {
    return nullptr;
  }
  SpatialObject new_entry;
  new_entry.box.reset();
  for (SpatialObject &current_entry : *right) {
    new_entry.box.adjust(current_entry.box);
  }
  new_entry.child_pointer = right;
  return parent->insert(new_entry);
}
///
template <size_t N, typename ElemType, size_t M, size_t m>
	void RTree<N, ElemType, M, m>::buscar(Rectangle<N>& rectangulo,
										 std::shared_ptr<Node> nodo_actual) {
	//primera vez
	if (!nodo_actual){
		nodo_actual = root_pointer_;
	} 
	
	for (const auto& ptr : *nodo_actual) {
		if (overlaps(ptr.box, rectangulo)) {
			if (nodo_actual->is_leaf())
				historial.push_back(ptr.identifier);
			else
				buscar(rectangulo, ptr.child_pointer);
		}
	}
	
}
#endif  // SOURCE_RTREE_HPP_
