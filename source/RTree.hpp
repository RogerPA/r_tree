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
		typedef SpatialObject* iterator;
		typedef const SpatialObject* const_iterator;

		iterator begin();
		iterator end();

		const_iterator begin() const;
		const_iterator end() const;

		SpatialObject& operator[](size_t index);
		SpatialObject operator[](size_t index) const;

		bool is_leaf();
		void pickSeeds(std::vector<SpatialObject>& lista, size_t* a, size_t* b);
		void pickNext(std::vector<SpatialObject>& list_entri, size_t* i_x, const Rectangle<N>& b_l, const Rectangle<N>& b_r);

		std::shared_ptr<Node> insert(const SpatialObject& new_entry);

		SpatialObject entry[M];
		size_t size = 0;
 };

  RTree();
  virtual ~RTree();
  size_t dimension() const;
  size_t size() const;
  bool empty() const;

  void insert(const Rectangle<N>& box, const ElemType& value);
  std::shared_ptr<Node> choose_leaf(const std::shared_ptr<Node>& current_node,
                                    const Rectangle<N>& box,
                                    const ElemType& value);

  std::shared_ptr<Node> choose_node(const std::shared_ptr<Node>& current_node,
                                    const Rectangle<N>& box,
                                    SpatialObject*& entry);

  std::shared_ptr<Node> adjust_tree(const std::shared_ptr<Node>& parent,
                                    const std::shared_ptr<Node>& left,
                                    const std::shared_ptr<Node>& right,
                                    SpatialObject* entry);

  std::vector<ElemType>& operator[](const Rectangle<N>& box);
  std::vector<ElemType>& at(const Rectangle<N>& box);
  const std::vector<ElemType>& at(const Rectangle<N>& box) const;
  void query(const Rectangle<N>& box, std::shared_ptr<Node> current = nullptr);
  // std::vector<ElemType> kNNValue(const Rectangle<N> &box, size_t k) const;

  // private:
  std::shared_ptr<Node> root_pointer_;
  size_t entry_count;
  std::vector<ElemType> query_result;
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
typename RTree<N, ElemType, M, m>::SpatialObject&
RTree<N, ElemType, M, m>::Node::operator[](size_t index) {
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
RTree<N, ElemType, M, m>::Node::insert(const SpatialObject& new_entry) {
	if (size < M) {
		entry[size++] = new_entry;
		return nullptr;
	}

	std::vector<SpatialObject> entries(M + 1);
	std::copy(begin(), end(), entries.begin());
	entries.at(M) = new_entry;
	this->size = 0;

	size_t first, second;
	pickSeeds(entries, &first, &second);

	this->insert(entries.at(first));
	auto new_node = std::make_shared<Node>();
	new_node->insert(entries.at(second));

	auto mbb_group1 = entries.at(first).box;
	auto mbb_group2 = entries.at(second).box;

	entries.erase(entries.begin() + first);
	entries.erase(entries.begin() + second - 1);

	float area_increase_g1, area_increase_g2;
	Rectangle<N> enlarged_mbb_g1, enlarged_mbb_g2;
	size_t idx;
	while (!entries.empty()) {
		if (this->size + entries.size() == m) {
			for (auto& e : entries) this->insert(e);
			entries.clear();
		} 
		else if (new_node->size + entries.size() == m) {
			for (auto& e : entries) new_node->insert(e);
			entries.clear();
		}
		else {
			pickNext(entries, &idx, mbb_group1, mbb_group2);

			enlarged_mbb_g1 = mbb_group1;
			enlarged_mbb_g1.adjust(entries.at(idx).box);
			area_increase_g1 = enlarged_mbb_g1.get_area() - mbb_group1.get_area();

			enlarged_mbb_g2 = mbb_group2;
			enlarged_mbb_g2.adjust(entries.at(idx).box);
			area_increase_g2 = enlarged_mbb_g2.get_area() - mbb_group2.get_area();

			if (area_increase_g1 < area_increase_g2 ||
			(area_increase_g1 == area_increase_g2 &&
			mbb_group1.get_area() < mbb_group2.get_area()) ||
			(area_increase_g1 == area_increase_g2 &&
			mbb_group1.get_area() == mbb_group2.get_area() &&
			size <= new_node->size)) {
				mbb_group1.adjust(entries.at(idx).box);
				this->insert(entries.at(idx));
			} 
			else {
				mbb_group2.adjust(entries.at(idx).box);
				new_node->insert(entries.at(idx));
			}

			entries.erase(entries.begin() + idx);
		}
	}

	return new_node;
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::Node::pickSeeds(std::vector<SpatialObject>& entries, size_t* first, size_t* second) {
	float max_waste = 0.f, area;
	Rectangle<N> enlargement_box;
	for (size_t i = 0; i < M; ++i) {
		for (size_t j = i + 1; j < M + 1; ++j) {
		enlargement_box = entries.at(i).box;
		enlargement_box.adjust(entries.at(j).box);
		area = enlargement_box.get_area() - entries.at(i).box.get_area() - entries.at(j).box.get_area();
		if (area > max_waste) {
			*first = i;
			*second = j;
		}
		}
	}
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::Node::pickNext(std::vector<SpatialObject>& entries, size_t* idx,
												const Rectangle<N>& mbb_group1, const Rectangle<N>& mbb_group2) {
	float max_diff = -1.f;
	Rectangle<N> enlarged_g1, enlarged_g2;
	float area_increase_g1, area_increase_g2, area_diff;
	for (size_t i = 0; i < entries.size(); ++i) {
		enlarged_g1 = mbb_group1;
		enlarged_g1.adjust(entries.at(i).box);
		area_increase_g1 = enlarged_g1.get_area() - mbb_group1.get_area();

		enlarged_g2 = mbb_group2;
		enlarged_g2.adjust(entries.at(i).box);
		area_increase_g2 = enlarged_g2.get_area() - mbb_group2.get_area();

		area_diff = std::abs(area_increase_g1 - area_increase_g2);

		if (area_diff > max_diff) {
			max_diff = area_diff;
			*idx = i;
		}
	}
}

/** R-Tree class implementation details */

template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::RTree() : root_pointer_(new Node), entry_count(0) {}

template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::~RTree() {}

template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::dimension() const {
	return N;
}

template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::size() const {
	return entry_count;
}

template <size_t N, typename ElemType, size_t M, size_t m>
bool RTree<N, ElemType, M, m>::empty() const {
	return size() == 0;
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::insert(const Rectangle<N>& box, const ElemType& value) {
	std::shared_ptr<Node> splitted_node = choose_leaf(root_pointer_, box, value);
	++entry_count;
	if (!splitted_node) {
		return;
	}

	auto new_root = std::make_shared<Node>();
	new_root->entry[0].child_pointer = root_pointer_;
	++new_root->size;
	adjust_tree(new_root, root_pointer_, splitted_node, &new_root->entry[0]);
	root_pointer_ = new_root;
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::choose_leaf(const std::shared_ptr<Node>& current_node, const Rectangle<N>& box, const ElemType& value) {
	if (!current_node->is_leaf()) {
		SpatialObject* entry;
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
RTree<N, ElemType, M, m>::choose_node(const std::shared_ptr<Node>& current_node, const Rectangle<N>& box, SpatialObject*& entry) {
	float minimum_area = (*current_node)[0].box.get_area();

	Rectangle<N> enlarged_box = (*current_node)[0].box;
	enlarged_box.adjust(box);
	float minimum_enlargement = enlarged_box.get_area() - minimum_area;

	float enlargement, area;
	std::shared_ptr<Node> node = (*current_node)[0].child_pointer;

	entry = &(*current_node)[0];
	for (SpatialObject& current_entry : *current_node) {
		area = current_entry.box.get_area();

		enlarged_box = current_entry.box;
		enlarged_box.adjust(box);
		enlargement = enlarged_box.get_area() - area;

		if (enlargement < minimum_enlargement || (enlargement == minimum_enlargement && area < minimum_area)) {
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
RTree<N, ElemType, M, m>::adjust_tree(const std::shared_ptr<Node>& parent,
                                      const std::shared_ptr<Node>& left,
                                      const std::shared_ptr<Node>& right,
                                      SpatialObject* entry) {
	entry->box.reset();
	for (SpatialObject current_entry : *left) {
		entry->box.adjust(current_entry.box);
	}
	if (!right) {
		return nullptr;
	}
	SpatialObject new_entry;
	new_entry.box.reset();
	for (SpatialObject& current_entry : *right) {
		new_entry.box.adjust(current_entry.box);
	}
	new_entry.child_pointer = right;
	return parent->insert(new_entry);
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::vector<ElemType>& RTree<N, ElemType, M, m>::operator[](const Rectangle<N>& box) {
	query_result.clear();
	query(box);
	return query_result;
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::vector<ElemType>& RTree<N, ElemType, M, m>::at(const Rectangle<N>& box) {
	query_result.clear();
	query(box);
	return query_result;
}

template <size_t N, typename ElemType, size_t M, size_t m>
const std::vector<ElemType>& RTree<N, ElemType, M, m>::at(const Rectangle<N>& box) const {
	query_result.clear();
	query(box);
	return query_result;
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::query(const Rectangle<N>& box, std::shared_ptr<Node> current) {
	if (!current) current = root_pointer_;

	for (const auto& so : *current) {
		if (overlaps(so.box, box)) {
			if (current->is_leaf())
				query_result.push_back(so.identifier);
			else
				query(box, so.child_pointer);
		}
	}
}

#endif  // SOURCE_RTREE_HPP_
