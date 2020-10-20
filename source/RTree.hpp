// Copyright

#ifndef SOURCE_RTREE_HPP_
#define SOURCE_RTREE_HPP_

#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

#include "Rectangle.hpp"

template <size_t N, typename ElemType, size_t M, size_t m = M / 2>//dimensiones,tipo de dato,max #, min #
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
    void pick_seeds(int& a,int& b);
    int pick_next(const int& i, const int& a, const int& b);

    SpatialObject entry[M];
    size_t size = 0;
  };

  RTree();
  virtual ~RTree();
  size_t dimension() const;
  size_t size() const;
  bool empty() const;

  void insert(const Rectangle<N> &box, const ElemType &value);

  std::shared_ptr<Node> choose_leaf(const std::shared_ptr<Node> &current_node,const Rectangle<N> &box,const ElemType &value);

  std::shared_ptr<Node> choose_node(const std::shared_ptr<Node> &current_node,const Rectangle<N> &box,SpatialObject *&entry);

  std::shared_ptr<Node> adjust_tree(const std::shared_ptr<Node> &parent,const std::shared_ptr<Node> &left,const std::shared_ptr<Node> &right,
                                    SpatialObject *entry);

  // TODO(ADE): Implement the details of all this functions
  std::vector<ElemType> &operator[](const Rectangle<N> &box);
  std::vector<ElemType> &at(const Rectangle<N> &box);
  const std::vector<ElemType> &at(const Rectangle<N> &box) const;
  // std::vector<ElemType> kNNValue(const Rectangle<N> &box, size_t k) const;

  // private:
  std::shared_ptr<Node> root_pointer_;
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
void RTree<N, ElemType, M, m>::Node::pick_seeds(int& a,int& b) {
    float most_pair = 0;float d = 0;
    Rectangle<N> aux = entry[0].box;

    for (int j = 0; j < M;j++) {
        for (int i = j + 1; i < M;i++) {
            aux.reset();
            aux.adjust(entry[j].box);
            aux.adjust(entry[i].box);

            d = (aux.get_area() - entry[j].box.get_area() - entry[i].box.get_area());

            if (d > most_pair) {
                most_pair = d;
                a = j; b = i;
            }

        }
    }
}

template<size_t N, typename ElemType, size_t M, size_t m>
int RTree<N, ElemType, M, m>::Node::pick_next(const int& i, const int& a,const int& b){
    float most = 0;float d = 0;float d2 = 0;

    for (int j = 0; j < N; j++) {
        d += get_enlargement(entry[a].box[j], entry[i].box[j]);
    }

    for (int j = 0; j < N; j++) {
        d2 += get_enlargement(entry[b].box[j], entry[i].box[j]);
    }

    if (d > d2) {
        return b;
    }
    else if (d2 > d) {
        return a;
    }

    Rectangle<N> aux = entry[a].box;
    Rectangle<N> aux2 = entry[b].box;

    aux.reset();
    aux.adjust(entry[i].box);
    aux.adjust(entry[a].box);
    d = aux.get_area() - entry[a].box.get_area();

    aux2.reset();
    aux2.adjust(entry[i].box);
    aux2.adjust(entry[b].box);
    d2 = aux2.get_area() - entry[b].box.get_area();

    if (d > d2) {
        return b;
    }

    return a;
}

// TODO(ADE): Split the entries and return a pointer to new node
// caused due to split.

template <size_t N, typename ElemType, size_t M, size_t m>
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::Node::insert(const SpatialObject &new_entry) {
  if (size + 1 <= M) {
    entry[size++] = new_entry;
    return nullptr;
  }
  
  int a = 0 ,b = 0;
  pick_seeds(a,b);
  std::shared_ptr<Node> right = std::make_shared<Node>();
  std::shared_ptr<Node> left = std::make_shared<Node>();

  for (int i = 0; i < M;i++) {
      if (i == a || i == b) {
          continue;
      }
      int c = pick_next(i, a, b);
      if (c == b) {
          right->insert(entry[i]);
      }
      left->insert(entry[i]);
  }

  //this = left;

  return right;
}


/** R-Tree class implementation details */

template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::RTree() : root_pointer_(new Node) {}

// TODO(ADE):
template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::~RTree() {
}

// TODO(ADE):
template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::dimension() const {
  return N;
}

// TODO(ADE):
template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::size() const {
  return size_t(0);
}

// TODO(ADE):
template <size_t N, typename ElemType, size_t M, size_t m>
bool RTree<N, ElemType, M, m>::empty() const {
    if (this->size() >= 1) {
        return false;
    }
    return true;
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::insert(const Rectangle<N> &box, const ElemType &value) {

  std::shared_ptr<Node> splitted_node = choose_leaf(root_pointer_, box, value);
  if (!splitted_node) {
    return;
  }
  // TODO(ADE): Last part of insert is missing i.e. when the root overflow
  // see R-tree gutman paper description.
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::choose_leaf(const std::shared_ptr<Node> &current_node,const Rectangle<N> &box,const ElemType &value) {

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
RTree<N, ElemType, M, m>::choose_node(const std::shared_ptr<Node> &current_node,const Rectangle<N> &box,SpatialObject *&entry) {
  //ELIGE A QUE NODO CONTINUAR EL CHOOSE LEAF
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

#endif  // SOURCE_RTREE_HPP_
