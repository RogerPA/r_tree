// Copyright

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

    std::shared_ptr<Node> insert(const SpatialObject &new_entry);
	
	void PickSeeds(std::vector<SpatialObject> &entries, size_t *first, size_t *second);
	void PickNext(std::vector<SpatialObject> &entries, size_t *iter, Rectangle<N> &mbb_1, Rectangle<N> &mbb_2);
	

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
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::Node::insert(const SpatialObject &new_entry) {
  if (size < M) {
    entry[size++] = new_entry;
    return nullptr;
  }
  // TODO(ADE): Split the entries and return a pointer to new node
  // caused due to split.
  std::vector<SpatialObject> entries(M+1);
  std::shared_ptr<Node> new_node = std::make_shared<Node>();
  for(size_t i = 0; i<M; ++i){
	  entries[i] = entry[i];
  }
  entries[M] = new_entry;
  this->size = 0;
  size_t first, second;
  //Pick first entry for each group
  PickSeeds(entries, &first, &second);
  //Assign each to a group
  insert(entries[first]);
  new_node->insert(entries[second]);
  if(first > second){
	  std::swap(first,second);
  }
  entries.erase(entries.begin() + first);
  entries.erase(entries.begin() + second - 1);
  Rectangle<N> d1 = (*this).box;
  Rectangle<N> d2 = (*new_node).box;
  Rectangle<N> enlargement_box_g1, enlargement_box_g2;
  float area_g1, area_g2;
  size_t iter;
  //Repeat from QS2 and QS3 if not empty
  while(!entries.empty()){
	  //Check if done
	  //If all entries have been assigned, stop
	  if(size + entries.size() == m){
		  for(auto &e : entries){
			  insert(e);
		  }
		  entries.clear();
		  continue;
	  }
	  //If one group has so few entries have been assigned in order, stop
	  else if(new_node->size + entries.size() == m){
		  for(auto &e : entries){
			  insert(e);
		  }
		  entries.clear();
		  continue;
	  }
	  //If not done select the entry to assign it
	  else{
		  //INvoke PickNext
		  PickNext(entries,iter,d1,d2);
		  //Rectangles that will have to be enlarged
		  enlargement_box_g1 = d1;
		  enlargement_box_g1.adjust(entries[iter].box);
		  area_g1 = enlargement_box_g1.get_area() - d1.get_area();
		  enlargement_box_g2 = d2;
		  enlargement_box_g2.adjust(entries[iter].box);
		  area_g2 = enlargement_box_g2.get_area() - d2.get_area();
		  //Resolve ties by adding the entry to the group with smallest area
		  if(area_g1 < area_g2 || 
			 (area_g1 == area_g2 && d1.get_area() < d2.get_area()) ||
			 (area_g1 == area_g2 && d1.get_area() == d2.get_area()) &&
			 size < new_node->size){
			  d1.adjust(entries[iter].box);
			  //Add it to the group
			  insert(entries[iter]);
		  }
		  else{
			  d2.adjust(entries[iter].box);
			  //Add it to the group
			  insert(entries[iter]);
		  }
		  entries.erase(entries.begin()+iter);
	  }
  }
  return new_node;
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::Node::PickSeeds(std::vector<SpatialObject> &entries, size_t *first, size_t *second){
	//Calculate inefficiency of grouping entries together
	Rectangle<N> enlargement_box;
	float max_area = 0.0f;
	for(size_t i = 0; i<M; ++i){
		for(size_t j = i+1; j<M; ++j){
			enlargement_box = entries[i].box;
			enlargement_box.adjust(entries[j].box);
			float d = enlargement_box.get_area() - entries[i].box.get_area() - entries[j].box.get_area();
			max_area = d;
			if(d > max_area){
				//most wasteful pair with the largest d(area)
				//max_area = d; 
				*first = i;
				*second = j;
			}
		}
	}
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::Node::PickNext(std::vector<SpatialObject> &entries, size_t *iter, Rectangle<N> &mbb_1, Rectangle<N> &mbb_2){
	//Determine the cost of putting each entry in each group
	Rectangle<N> rectangle_group1, rectangle_group2;
	float d1, d2; //area increase required in the covering rectangle of each group
	float max_difference_area = 0.0f;
	for(size_t i = 0; i<entries.size(); ++i){
		rectangle_group1 = mbb_1;
		rectangle_group1.adjust(entries[i].box);
		d1 = rectangle_group1.get_area() - mbb_1.get_area();
		rectangle_group2 = mbb_2;
		rectangle_group2.adjust(entries[j].box);
		d2 = rectangle_group2.get_area() - mbb_2.get_area();
		//Find entry with greatest preference for one group
		float d = (d1 - d2);
		if(d > max_difference_area){
			max_difference_area = d;
			*iter = i;
		}
	}
}

/** R-Tree class implementation details */

template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::RTree() : root_pointer_(new Node) {}

// TODO(ADE):
template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::~RTree() {}

// TODO(ADE):
template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::dimension() const {
  return size_t(0);
}

// TODO(ADE):
template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::size() const {
  return size_t(0);
}

// TODO(ADE):
template <size_t N, typename ElemType, size_t M, size_t m>
bool RTree<N, ElemType, M, m>::empty() const {
  return false;
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
