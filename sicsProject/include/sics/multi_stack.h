#ifndef SICS_MULTI_STACK_H_
#define SICS_MULTI_STACK_H_

#include <iterator>
#include <vector>

namespace sics {

template <typename T>
class multi_stack {
 public:
  using iterator = typename std::vector<T>::reverse_iterator;
  using size_type = typename std::vector<T>::size_type;

 private:
  std::vector<T> data;
  std::vector<iterator> level_its;
  iterator it;
  typename std::vector<iterator>::reverse_iterator level_it;

 public:
  multi_stack(size_type n, size_type m)
      : data(n),
        level_its(m),
        it{std::rend(data)},
        level_it{std::rend(level_its)} {
  }

  multi_stack(multi_stack const &) = delete;

  bool empty() const {
    return it == std::rend(data);
  }

  auto size() const {
    return std::distance(it, std::rend(data));
  }

  T const & top() const {
    return *it;
  }

  void pop() {
    ++it;
  }

  void push(T item) {
    --it;
    *it = item;
  }

  bool level_empty() const {
    return it == *level_it;
  }

  auto level_size() const {
    return std::distance(it, *level_it);
  }

  void push_level() {
    --level_it;
    *level_it = it;
  }

  void pop_level() {
    ++level_it;
  }

};

}  // namespace sics

#endif  // SICS_MULTI_STACK_H_
