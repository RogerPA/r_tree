// Copyright

#ifndef SOURCE_RTREE_HPP_
#define SOURCE_RTREE_HPP_

#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>
#include <utility>
#include <math.h>

#include "Rectangle.hpp"
//DIMENSION - TIPO DE ELEMENTO - ESPACIOS POR NODO COMO MAXIMO - ESPACIOS POR NODO COMO MINIMO
template <size_t N, typename ElemType, size_t M, size_t m = M / 2>
class RTree {
public:
  struct Node;

  struct SpatialObject {
    Rectangle<N> box;
    ElemType identifier;
    std::shared_ptr<Node> child_pointer;

    SpatialObject& operator=(const SpatialObject& spatial) {
      box = spatial.box;
      identifier = spatial.identifier;
      child_pointer = spatial.child_pointer;
      vi2split = false;
    }

    bool vi2split = false;
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
    std::pair<std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>, std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>>
                                         split_node(const SpatialObject &new_entry);
    std::pair<size_t, size_t> pick_seeds();
    std::pair<size_t, bool> pick_next(std::vector<SpatialObject> &freeEntries,
                                         std::pair<Rectangle<N>, Rectangle<N>> &rect1_2,
                                         std::shared_ptr<Node> g1, std::shared_ptr<Node> g2);

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
  // std::vector<ElemType> kNNValue(const Rectangle<N> &box, size_t k) const;

//private:
  std::shared_ptr<Node> root_pointer_;
  size_t num_of_elements;
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
std::pair<size_t,size_t> RTree<N, ElemType, M, m>::Node::pick_seeds() {
  float d = 0.0;
  std::pair<size_t, size_t> seeds; seeds.first = seeds.second = size_t(0);
  for (size_t sp1 = 0; sp1 < M; sp1++) {
    for (size_t sp2 = sp1 + 1; sp2 < M; sp2++) {
      Rectangle<N> j = entry[sp1].box;
      j.adjust(entry[sp2].box);
      float temp = d;
      d = max(abs(j.get_area() - entry[sp1].box.get_area() - entry[sp2].box.get_area()), d);
      if (temp != d)
        seeds = std::make_pair(sp1, sp2);
    }
  }
  return seeds;
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::pair<size_t, bool> RTree<N, ElemType, M, m>::Node::pick_next(std::vector<SpatialObject> &freeEntries,
                                                                 std::pair<Rectangle<N>,Rectangle<N>> &rect1_2,
                                                                 std::shared_ptr<Node> g1, std::shared_ptr<Node> g2) {
  float d1 = 0.0, d2 = 0.0, difmax = 0.0;
  std::pair<size_t, bool> newEntry = std::make_pair(size_t(0),0);
  for (size_t i = 0; i < freeEntries.size(); i++) {
    rect1_2.first.adjust(freeEntries[i].box);
    rect1_2.second.adjust(freeEntries[i].box);
    d1 = abs(rect1_2.first.get_area() - freeEntries[i].box.get_area());
    d2 = abs(rect1_2.second.get_area() - freeEntries[i].box.get_area());
    float temp = difmax;
    difmax = max(abs(d1 - d2),difmax);
    if (temp != difmax) {
      newEntry.first = i;
      newEntry.second = (d1 > d2);//0-> go to g1, 1-> go to g2
    }
  }
  return newEntry;
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::pair<std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>, std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>>
RTree<N, ElemType, M, m>::Node::split_node(const SpatialObject &new_entry) {
  std::pair<size_t,size_t> seeds = pick_seeds();

  std::shared_ptr<Node> group1(new Node), group2(new Node);
  group1->entry[group1->size++] = entry[seeds.first]; entry[seeds.first].vi2split = true;
  group2->entry[group2->size++] = entry[seeds.second]; entry[seeds.second].vi2split = true;

  std::vector<SpatialObject> lastgroup; lastgroup.resize(M - 1);//M(entries who are already in the node) + 1(new entry) - 2(seeds)
  size_t i = size_t(0);
  for (SpatialObject &sp_obj : entry) {
    if (!sp_obj.vi2split) {
      sp_obj.vi2split = false;
      lastgroup[i++] = sp_obj;
    }
  }
  lastgroup[i] = new_entry;

  while (!lastgroup.empty()) {
    std::pair<size_t, bool> an_entry = pick_next(lastgroup, group1, group2);
    if (an_entry.second) {
      group1->entry[group1->size++] = lastgroup[an_entry.first];
    }
    else {
      group2->entry[group2->size++] = lastgroup[an_entry.first];
    }
    lastgroup.erase(lastgroup.begin() + an_entry.first);
  }
  if (group1->size < m) {
    SpatialObject voidObj;
    group1->entry[group1->size++] = (*group2)[group2->size];
    (*group2)[group2->size--] = voidObj;
  }
  else if (group2->size < m) {
    SpatialObject voidObj;
    group2->entry[group2->size++] = (*group1)[group1->size];
    (*group1)[group1->size--] = voidObj;
  }
  return std::make_pair(group1,group2);
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::Node::insert(const SpatialObject &new_entry) {
  if (size < M) {
    entry[size++] = new_entry;
    return nullptr;
  }
  // TODO(ADE): Split the entries and return a pointer to new node
  // caused due to split.
  split_node(new_entry);

  return nullptr;
}

/** R-Tree class implementation details */

// TODO(ADE):(COMPLETE)
template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::RTree() : root_pointer_(new Node) { num_of_elements = size_t(0); }

// TODO(ADE):(COMPLETE ????)
template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::~RTree() { }

// TODO(ADE):(COMPLETE)
template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::dimension() const {
  return N;
}

// TODO(ADE):(COMPLETE)
template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::size() const {
  return num_of_elements;
}

// TODO(ADE):(COMPLETE)
template <size_t N, typename ElemType, size_t M, size_t m>
bool RTree<N, ElemType, M, m>::empty() const {
  return !root_pointer_->size;
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::insert(const Rectangle<N> &box,
                                      const ElemType &value) {
  std::shared_ptr<Node> splitted_node = choose_leaf(root_pointer_, box, value);
  if (!splitted_node) {
    return;
  }
  // TODO(ADE): Last part of insert is missing i.e. when the root overflow
  // see R-tree gutman paper description.
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
  //std::cout << "Current minimum" << entry->box << std::endl;
  for (SpatialObject &current_entry : *current_node) {
    //std::cout << "Comparing with " << current_entry.box << std::endl;
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