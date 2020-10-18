// Copyright

#ifndef SOURCE_RTREE_HPP_
#define SOURCE_RTREE_HPP_

#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>
#include <queue>

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

		void elegir_semilla(std::vector<SpatialObject> &lista, int &a, int &b);
		int elegir_siguiente(std::vector<SpatialObject> &lista_, Rectangle<N> &b_left, Rectangle<N> &b_right);
    
		SpatialObject entry[M];
		size_t size = 0;
	};

	RTree();
	virtual ~RTree();
	size_t dimension() const;
	size_t size() const;
	bool empty() const;
	void print();

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

	std::vector<ElemType> operator[](const Rectangle<N> &box);
	std::vector<ElemType> at(const Rectangle<N> &box);
	const std::vector<ElemType> at(const Rectangle<N> &box) const;
	// std::vector<ElemType> kNNValue(const Rectangle<N> &box, size_t k) const;

	void buscar(const std::shared_ptr<Node> &cur_node,
				const Rectangle<N> &box,
				std::vector<ElemType> &result);

	//private:
	std::shared_ptr<Node> root_pointer_;
	size_t entries;
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
  
	 // Split the entries and return a pointer to new node caused due to split.
	std::shared_ptr<Node> new_node = std::make_shared<Node>();
	std::vector<SpatialObject> entry_list(M + 1);
	for (size_t i = 0; i < M; ++i) {
		entry_list[i] = entry[i];
	}
	entry_list[M] = new_entry;
	this->size = 0;
  
	int left, right;
	elegir_semilla(entry_list, left, right);
	this->insert(entry_list[left]);
	new_node->insert(entry_list[right]);
	if(left > right) std::swap(left, right);
	entry_list.erase(entry_list.begin() + left);
	entry_list.erase(entry_list.begin() + right - 1);
  
	Rectangle<N> bb_left = (*this)[0].box;
	Rectangle<N> bb_right = (*new_node)[0].box;
	while (!entry_list.empty()) {
		if (this->size + entry_list.size() == m) {
			for (SpatialObject &obj : entry_list) {
				this->insert(obj);
			}
			entry_list.clear();
			continue;
		}
		if (new_node->size + entry_list.size() == m) {
			for (SpatialObject &obj : entry_list) {
				new_node->insert(obj);
			}
		entry_list.clear();
		continue;
		}
		int next = elegir_siguiente(entry_list, bb_left, bb_right);
		Rectangle<N> enlargement_left = entry_list[next].box;
		Rectangle<N> enlargement_right = entry_list[next].box;
		enlargement_left.adjust(bb_left);
		enlargement_right.adjust(bb_right);
		float area_left = enlargement_left.get_area() - entry_list[next].box.get_area();
		float area_right = enlargement_right.get_area() - entry_list[next].box.get_area();
		if (area_left < area_right ||
        (area_left == area_right && bb_left.get_area() < bb_right.get_area()) ||
        (area_left == area_right && bb_left.get_area() == bb_right.get_area() &&
         this->size < new_node->size)) {
			bb_left.adjust(entry_list[next].box);
			this->insert(entry_list[next]);
		}
		else {
			bb_right.adjust(entry_list[next].box);
			new_node->insert(entry_list[next]);
		}
		entry_list.erase(entry_list.begin() + next);
	}
	return new_node;
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::Node::elegir_semilla(std::vector<SpatialObject> &list_, int &a, int &b) {
	a = 0, b = 1;
	Rectangle<N> enlargement_box = list_[a].box;
	enlargement_box.adjust(list_[b].box);
	float max_area = enlargement_box.get_area();
	max_area -= (list_[a].box.get_area() + list_[b].box.get_area());
	float area;
	for (int i = 0; i < M; ++i) {
		for (int j = i+1; j <= M; ++j) {
			enlargement_box = list_[i].box;
			enlargement_box.adjust(list_[j].box);
			area = enlargement_box.get_area() - list_[i].box.get_area() - list_[j].box.get_area();
			if (area > max_area) {
				max_area = area;
				a = i, b = j;
			}
		}
	}
}

template <size_t N, typename ElemType, size_t M, size_t m>
int RTree<N, ElemType, M, m>::Node::elegir_siguiente(std::vector<SpatialObject> &list_, Rectangle<N> &bb_left, Rectangle<N> &bb_right) {
	int ind = 0;
	Rectangle<N> enlargement_left = list_[0].box;
	Rectangle<N> enlargement_right = list_[0].box;
	enlargement_left.adjust(bb_left);
	enlargement_right.adjust(bb_right);
	float area_left = enlargement_left.get_area() - list_[0].box.get_area();
	float area_right = enlargement_right.get_area() - list_[0].box.get_area();
	float max_d = abs(area_left - area_right);
	float area_d;
	for (int i = 1; i < list_.size(); ++i) {
		enlargement_left = list_[i].box;
		enlargement_right = list_[i].box;
		enlargement_left.adjust(bb_left);
		enlargement_right.adjust(bb_right);
		area_left = enlargement_left.get_area() - list_[i].box.get_area();
		area_right = enlargement_right.get_area() - list_[i].box.get_area();
		area_d = abs(area_left - area_right);
		if (area_d > max_d) {
		  max_d = area_d;
		  ind = i;
		}
	}
	return ind;
}

/** R-Tree class implementation details */

template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::RTree() : root_pointer_(std::make_shared<Node>()), entries(0) {}

template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::~RTree() { root_pointer_.reset(); }

template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::dimension() const {
	return N;
}

template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::size() const {
	return entries;
}

template <size_t N, typename ElemType, size_t M, size_t m>
bool RTree<N, ElemType, M, m>::empty() const {
	return !entries;
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::insert(const Rectangle<N> &box, const ElemType &value) {
	std::shared_ptr<Node> splitted_node = choose_leaf(root_pointer_, box, value);
	if (!splitted_node) {
		return;
	}
	// Split the root !
	std::shared_ptr<Node> new_root = std::make_shared<Node>();
	(*new_root)[0].child_pointer = root_pointer_;
	(new_root->size)++;
	adjust_tree(new_root, root_pointer_, splitted_node, &(*new_root)[0]);
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
	new_entry.child_pointer = nullptr;
	entries++;
	return current_node->insert(new_entry);
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::choose_node(const std::shared_ptr<Node> &current_node,
                                      const Rectangle<N> &box,
                                      SpatialObject *&entry) {
	Rectangle<N> enlarged_box = (*current_node)[0].box;
	enlarged_box.adjust(box);
	float minimum_area = (*current_node)[0].box.get_area();
	float minimum_enlargement = enlarged_box.get_area() - minimum_area;
	std::shared_ptr<Node> node = (*current_node)[0].child_pointer;
	float enlargement, area;

	entry = &(*current_node)[0];
	for (SpatialObject &current_entry : *current_node) {
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

template <size_t N, typename ElemType, size_t M, size_t m>
std::vector<ElemType> RTree<N, ElemType, M, m>::operator[](const Rectangle<N> &box) {
	std::vector<ElemType> result;
	search(root_pointer_, box, result);
	return result;
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::vector<ElemType> RTree<N, ElemType, M, m>::at(const Rectangle<N> &box){
	std::vector<ElemType> result;
	search(root_pointer_, box, result);
	return result;
}

template <size_t N, typename ElemType, size_t M, size_t m>
const std::vector<ElemType> RTree<N, ElemType, M, m>::at(const Rectangle<N> &box) const {
	std::vector<ElemType> result;
	buscr(root_pointer_, box, result);
	return result;
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::search(const std::shared_ptr<Node> &cur_node,
                                      const Rectangle<N> &box,
                                      std::vector<ElemType> &result){
	for (SpatialObject &obj : *cur_node) {
		if (overlaps(obj.box, box)) {
			if (cur_node->is_leaf()) {
				result.push_back(obj.identifier);
			}
		else {
			buscar(obj.child_pointer, box, result);
		}
		}
	}
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::print() {
	std::cout<<"Printing the tree ...";
	int lastLvl = -1;
	std::queue<std::pair<std::shared_ptr<Node>, int> > Q;
	Q.push(make_pair(root_pointer_, 0));
	while (!Q.empty()) {
		std::shared_ptr<Node> current = Q.front().first;
		int lvl = Q.front().second;
		Q.pop();
		if (lvl != lastLvl) {
			std::cout<<std::endl;
			lastLvl = lvl;
		}
		std::cout<<"{ | ";
		for (SpatialObject &obj : *current) {
			std::cout<<(obj.identifier.size() ? obj.identifier : "IN");
			std::cout<<":("<<obj.box[0].begin()<<","<<obj.box[0].end()<<") - ";
			std::cout<<"("<<obj.box[1].begin()<<","<<obj.box[1].end()<<") | ";
			if (!current->is_leaf()) {
				Q.push(make_pair(obj.child_pointer, lvl+1));
			}
		}
		std::cout<<"} ";    
	}
	std::cout<<"\n";
}

#endif  // SOURCE_RTREE_HPP_
