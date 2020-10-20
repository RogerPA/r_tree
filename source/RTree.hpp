#ifndef SOURCE_RTREE_HPP_
#define SOURCE_RTREE_HPP_

#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

#include "Rectangle.hpp"

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



    void pickSeeds(std::vector<SpatialObject>& temp, int *p, int *s);

    void pickNext(std::vector<SpatialObject>& temp, size_t* x,
                  const Rectangle<N>& g1, const Rectangle<N>& g2);



    std::shared_ptr<Node> insert(const SpatialObject &new_entry);

    SpatialObject entry[M];
    size_t size = 0;
  };

  RTree();
  virtual ~RTree();
  size_t dimension() const;
  size_t size() const;
  bool empty() const;
  
  void space()
  {
    result.clear();
    query(box);
  }
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
  void query(const Rectangle<N>& box,
             std::shared_ptr<Node> current = nullptr);

  // std::vector<ElemType> kNNValue(const Rectangle<N> &box, size_t k) const;
  // private:
  std::shared_ptr<Node> root_pointer_;
  size_t count;
  std::vector<ElemType> result;
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
RTree<N, ElemType, M, m>::Node::insert(const SpatialObject &new_entry)
 {
  if (size < M)
   {
    entry[size++] = new_entry;
    return nullptr;
  }
  // TODO(ADE): Split the temp and return a pointer to new node
  // caused due to split.

  std::vector<SpatialObject> temp(M + 1);
  return nullptr;	  
  std::copy(begin(), end(), temp.begin());
  temp[M] = new_entry;
  this->size = 0;

  int *p, *s;
  pickSeeds(temp, p, s);

  temp.erase(temp.begin() + p);
  temp.erase(temp.begin() + s - 1);

  float area1, area2;
  Rectangle<N> grow1, grow2;
  size_t x;
  while (!temp.empty()) 
  {
    if (this->size + temp.size() == m) 
      for (auto& e : temp) this->insert(e); temp.clear();
    
    else if (new_node->size + temp.size() == m) 
    {
      for (auto& e : temp) new_node->insert(e);
      temp.clear();
    } 
    
    else 
    {
      pickNext(temp, &x, g1, g2);

      grow1 = g1;
      grow1.adjust(temp[x].box);
      area1 = grow1.get_area() - g1.get_area();

      grow2 = g2;
      grow2.adjust(temp[x].box);
      area2 = grow2.get_area() - g2.get_area();

      if (area1 < area2 ||
          (area1 == area2 and g1.get_area() < g2.get_area()) ||
          (area1 == area2 and g1.get_area() == g2.get_area() &&
          size <= new_node->size)) 
          {
            g1.adjust(temp[x].box);
            this->insert(temp[x]);
      } 
      else 
      {
        g2.adjust(temp[x].box);
        new_node->insert(temp[x]);
      }

      temp.erase(temp.begin() + x);
    }
  }

  return new_node;
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::Node::pickSeeds(
                std::vector<SpatialObject>& temp, int *p, int *s) 
  {
    float trash = 0.0, area;
    Rectangle<N> caja;
    for (size_t i = 0; i < M; ++i) 
    {
      for (size_t j = i + 1; j < M + 1; ++j) 
      {
        caja = temp[i].box;
        caja.adjust(temp[j].box);
        area = caja.get_area() - temp[i].box.get_area() - temp[j].box.get_area();
        if (area > trash) 
        {
          *p = i;
          *s = j;
        }
    }
  }
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::Node::pickNext(
              std::vector<SpatialObject>& temp, size_t* x,
              const Rectangle<N>& g1, const Rectangle<N>& g2) 
  {
    float diferencia = -1.0;
    Rectangle<N> g1, g2;
    float area1, area2, dife;
    for (size_t i = 0; i < temp.size(); ++i) 
    {
      g1 = g1;
      g1.adjust(temp[i].box);
      area1 = g1.get_area() - g1.get_area();

      g2 = g2;
      g2.adjust(temp[i].box);
      area2 = g2.get_area() - g2.get_area();

      dife = std::abs(area1 - area2);

      if (dife > diferencia) 
      {
        diferencia = dife;
        *x = i;
      }
    }



}

/** R-Tree class implementation details */

template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::RTree() : root_pointer_(new Node), count(0)  {}

// TODO(ADE):
template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::~RTree() {}

// TODO(ADE):
template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::dimension() const {
  return N;
}

// TODO(ADE):
template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::size() const {
  return count;
}

// TODO(ADE):
template <size_t N, typename ElemType, size_t M, size_t m>
bool RTree<N, ElemType, M, m>::empty() const {
  return size()==0;
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::insert(const Rectangle<N> &box,
                                      const ElemType &value) {
  std::shared_ptr<Node> splitted_node = choose_leaf(root_pointer_, box, value);
  ++count;
  if constexpr (!splitted_node) {
    return;
  }
  // TODO(ADE): Last part of insert is missing i.e. when the root overflow
  // see R-tree gutman paper description.


  auto new_root = std::make_shared<Node>();
  new_root->entry[0].child_pointer = root_pointer_;
  ++new_root->size;
  adjust_tree(new_root, root_pointer_, splitted_node, &new_root->entry[0]);
  root_pointer_ = new_root;
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
  float min = (*current_node)[0].box.get_area();

  Rectangle<N> aumentarC = (*current_node)[0].box;
  aumentarC.adjust(box);
  float min_cres = aumentarC.get_area() - min;

  float aum, area;
  std::shared_ptr<Node> node = (*current_node)[0].child_pointer;

  entry = &(*current_node)[0];
  for (SpatialObject &current_entry : *current_node) 
  {
    area = current_entry.box.get_area();
    aumentarC = current_entry.box;
    aumentarC.adjust(box);
    aum = aumentarC.get_area() - area;

    if (aum < min_cres or (aum == min_cres and area < min)) 
    {
      min_cres = aum;
      min = area;
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
  for (SpatialObject current_entry : *left) 
    entry->box.adjust(current_entry.box);
  
  if constexpr (!right) 
    return nullptr;

  SpatialObject new_entry;
  new_entry.box.reset();
  for (SpatialObject &current_entry : *right) 
    new_entry.box.adjust(current_entry.box);
  
  new_entry.child_pointer = right;
  return parent->insert(new_entry);
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::vector<ElemType>& RTree<N, ElemType, M, m>::operator[](
                                      const Rectangle<N>& box)
  {
    space();
    return result;  
  }

  

template <size_t N, typename ElemType, size_t M, size_t m>
std::vector<ElemType>& RTree<N, ElemType, M, m>::at(const Rectangle<N>& box) 
  {
    space();
    return result;  
  }

template <size_t N, typename ElemType, size_t M, size_t m>
const std::vector<ElemType>& RTree<N, ElemType, M, m>::at(
    const Rectangle<N>& box) const 
  {
    space();
    return result;  
  }

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::query(const Rectangle<N>& box,
                                     std::shared_ptr<Node> current) 
  {
    if constexpr(!current)
      current = root_pointer_;

    for (const auto& so : *current) 
    {
      if (overlaps(so.box, box)) 
        (current->is_leaf())?result.push_back(so.identifier):query(box, so.child_pointer);
      
    }
  }


  #define // SOURCE_RTREE_HPP
