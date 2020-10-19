// Copyright

#ifndef SOURCE_RTREE_HPP_
#define SOURCE_RTREE_HPP_

#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

#include "Rectangle.hpp"

// N es numero de dimensiones 
// ElemType Tidpo de dato a insertar
// M maximo numero que peude tener mi nodo
// m es el minimo que puede tener mi nodo
template <size_t N, typename ElemType, size_t M, size_t m = M / 2>
class RTree {
 public:
  struct Node;

// Box es el delimitador del Obejto
// tiene Identificador -> dato o registro
// child_pointer es un Puntero a otro Nodo 
// Aclaración:Si el identificador existe significa que es parte del nodo hoja entonces el child_pointer no apunta a nada
  struct SpatialObject {
    Rectangle<N> box;
    ElemType identifier;
    std::shared_ptr<Node> child_pointer;
  };

// Nodo va tener varios Obejtos Espaciales (SpatialObject)
// Aclaración nuestro nodo no apunta a otro nodo como las listas enlazadas sino que el obejto espacial es que apunta a otros nodos
// size -> cantidad de obejtos espaciales que tiene el nodo ,este atributo nos ayuda a determinar si eres un nodo hoja 
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

  // -------------
    //Intento 
    void Pick_Seeds(std::vector<SpatialObject> &ListasEspaciales, int &e1, int &e2 );

    void Pick_Next(std::vector<SpatialObject> &ListasEspaciales,Rectangle<N> &rect1, Rectangle<N> &rect2, int &e, int &ng);

  // -----------------
    std::shared_ptr<Node> insert(const SpatialObject &new_entry);

    SpatialObject entry[M];
    size_t size = 0;
  };


//Metodos del R-tree 
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
  
  // Buscar el vecino mas cercano
  // std::vector<ElemType> kNNValue(const Rectangle<N> &box, size_t k) const;

  // private:
  std::shared_ptr<Node> root_pointer_;

  size_t cantidad_elementos = 0;
};

// Metodos del Nodo
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

// Si tiene una netrada y una de esas entradas tiene puntero no es un hoja
template <size_t N, typename ElemType, size_t M, size_t m>
bool RTree<N, ElemType, M, m>::Node::is_leaf() {
  if (size && entry[0].child_pointer) {
    return false;
  }
  return true;
}

// Pick Begin
template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::Node::Pick_Seeds(std::vector<SpatialObject> &ListasEspaciales, int &e1, int &e2) {
  float d = 0;
  float temp = 0;

  for (int i = 0; i < M + 1; i++) {
    for (int j = i + 1; j < M + 1; j++) {
          Rectangle<N> jj = ListasEspaciales[i].box;
          jj.adjust(ListasEspaciales[j].box);
          temp = abs(jj.get_area() - ListasEspaciales[i].box.get_area() - ListasEspaciales[j].box.get_area());
          if (d < temp)
          {
            seedJ = i;
            seedJ = j;
            d = temp;
          }
    }
  }
}



template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::Node::Pick_Next(std::vector<SpatialObject> &ListasEspaciales, Rectangle<N> &rect1, Rectangle<N> &rect2, int &e, int &ng) {
  
  float d1;
  float d2;
  float ms;
  float temp = 0;

  for (int i = 0; i < Entries.size(); i++) {
      rect1.adjust(Entries[i].box);
      rect2.adjust(Entries[i].box);  
      
      d1 = abs(rect1.get_area() - Entries[i].box.get_area());
      d2 = abs(rect2.get_area() - Entries[i].box.get_area());

      temp = abs(d1 - d2);

      if(ms > temp){
        ms = temp;
        e = i;  
        if(d1 > d2)
          ng = 0; 
        else
         ng = 1;  
      }
  }
}
// Pick Fin




// Hacer (Falta) ,Aqui usar pickseeds y picknext
// El Algoritmo Cuadratico
template <size_t N, typename ElemType, size_t M, size_t m>
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::Node::insert(const SpatialObject &new_entry) {
  if (size < M) {
    entry[size++] = new_entry;
    return nullptr;
  }
  // TODO(ADE): Split the entries and return a pointer to new node
  // caused due to split.

  // Split y retornar el nuevo nodo

  //PickSeeds
  //PickNext

  // if (is necesario ahcer ajustar)
    // adjust_tree(ListaEntradas);

  //return nuevo nodo 
  return nullptr;
}



/** R-Tree class implementation details */
template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::RTree() : root_pointer_(new Node) {}

// Hacer -> (Falta) Implementar Destructor
template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::~RTree() {}

// Hacer -> Retornar la dimension
template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::dimension() const {
  return N;
}

// Hacer -> (falta) Retornar el tamaño
template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::size() const {
  return cantidad_elementos;
}

// Hacer -> Retorna verdadero si esta vacio
template <size_t N, typename ElemType, size_t M, size_t m>
bool RTree<N, ElemType, M, m>::empty() const {
   return !root_pointer_->size;
}

// Recibimos Box que tendra N dimensiones (intervalos)
// values en el ejemplo seria string
template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::insert(const Rectangle<N> &box,
                                      const ElemType &value) {

  std::shared_ptr<Node> splitted_node = choose_leaf(root_pointer_, box, value);
  if (!splitted_node) {
    return;
  }

  // TODO(ADE): Last part of insert is missing i.e. when the root overflow
  // see R-tree gutman paper description.

  std::shared_ptr<Node> new_root = std::make_shared<Node>();
  std::shared_ptr<Node> raiz_nueva = std::make_shared<Node>();
  (*raiz_nueva)[0].child_pointer = root_pointer_;
  (raiz_nueva->size)++;
  
  adjust_tree(raiz_nueva, root_pointer_, splitted_node, &(*raiz_nueva)[0]);
  root_pointer_ = raiz_nueva;
}

// Retorna el nodo hoja donde pondremos la entrada E
// Current_node inicia como raiz, y este ira actualizandose
template <size_t N, typename ElemType, size_t M, size_t m>
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::choose_leaf(const std::shared_ptr<Node> &current_node,
                                      const Rectangle<N> &box,
                                      const ElemType &value) {
  // Si no eres una hoja hay que seguir recorriendo el árbol
  if (!current_node->is_leaf()) {
    SpatialObject *entry;
    std::shared_ptr<Node> next_node = choose_node(current_node, box, entry);
    std::shared_ptr<Node> splitted_node = choose_leaf(next_node, box, value);
    return adjust_tree(current_node, next_node, splitted_node, entry);
  }
  // Eres nodo hoja y debemos insertar
  SpatialObject new_entry;
  new_entry.box = box;
  new_entry.identifier = value;
  return current_node->insert(new_entry);
}

// Retorna un puntero al nodo donde debo insertar el objeto criterio de estiramiento minimo
// entry ->(Objeto espacial) indice a que objeto espacial debo ir debido a que el objeto espacial solo es un puntero a otro nodo mas no al obejto espacial especifico ,
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


// N-> left 
// NN -> right
// Ajustar al padre
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

  // Retorna nulo si no hubo split 
  return parent->insert(new_entry);
}

#endif  // SOURCE_RTREE_HPP_
