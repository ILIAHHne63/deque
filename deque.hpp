#include <algorithm>
#include <exception>
#include <iostream>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <vector>

const int64_t kSize = 100;

template <typename T, typename Allocator = std::allocator<T>>
class Deque {
 public:
  using allocator_type = Allocator;
  using value_type = typename allocator_type::value_type;
  using alloc_traits = std::allocator_traits<Allocator>;

  template <bool IsConst>
  class common_iterator;

  template <bool IsConst>
  class common_reverse_iterator;

  template <bool IsConst>
  class common_iterator {
   public:
    using value_type = typename std::conditional<IsConst, const T, T>::type;
    using pointer = typename std::conditional<IsConst, const T*, T*>::type;
    using iterator_category = std::random_access_iterator_tag;
    using reference = typename std::conditional<IsConst, const T&, T&>::type;
    using difference_type = std::ptrdiff_t;

   private:
    T* pointer_;
    Deque<T, Allocator>* deque_;
    int64_t digit_pointer_;
    int64_t pointer_array_;

   public:
    common_iterator(Deque<T, Allocator>* deque) : deque_(deque) {
      if ((*deque_).empty()) {
        digit_pointer_ = -1;
        pointer_array_ = -1;
      }
    }

    common_iterator(T* pointer, Deque<T, Allocator>* deque,
                    int64_t digit_pointer, int64_t pointer_array);

    common_iterator(const common_iterator& other) { *this = other; }

    reference operator*() const { return *pointer_; }

    reference operator*() { return *pointer_; }

    int64_t update() { return pointer_array_ + ((*deque_).array_.size() / 4); }

    int64_t update() const {
      return pointer_array_ + ((*deque_).array_.size() / 4);
    }

    int64_t place();

    int64_t place() const;

    bool end_pointer();

    bool set() {
      return (digit_pointer_ == -1 && pointer_array_ == -1) || end_pointer();
    }

    void operators(int64_t& new_digit, common_iterator& iter);

    void operators(int64_t& new_digit, common_iterator& iter) const;

    common_iterator& operator++();

    common_iterator& operator--();

    void function_update();

    common_iterator operator--(int);

    common_iterator& operator=(const common_iterator& other);

    common_iterator operator++(int);

    common_iterator& operator+=(int64_t digit);

    common_iterator& operator-=(int64_t digit);

    common_iterator operator+(int64_t digit) const;

    common_iterator operator-(int64_t digit) const;

    difference_type operator-(const common_iterator& other) const;

    bool operator>=(const common_iterator& compare_iter) const {
      return place() >= compare_iter.place();
    }

    bool operator<=(const common_iterator& compare_iter) const {
      return place() <= compare_iter.place();
    }

    bool operator>(const common_iterator& compare_iter) const {
      return place() > compare_iter.place();
    }

    bool operator<(const common_iterator& compare_iter) const {
      return place() < compare_iter.place();
    }

    bool operator==(const common_iterator& compare_iter) const {
      return place() == compare_iter.place();
    }

    bool operator!=(const common_iterator& compare_iter) const {
      return place() != compare_iter.place();
    }

    int& index() { return place(); }

    typename std::conditional<IsConst, const T*, T*>::type operator->() const {
      return pointer_;
    }

    typename std::conditional<IsConst, const T*, T*>::type operator->() {
      return pointer_;
    }
  };

  using iterator = common_iterator<false>;
  using const_iterator = common_iterator<true>;
  using reverse_iterator = std::reverse_iterator<common_iterator<false>>;
  using const_reverse_iterator = std::reverse_iterator<common_iterator<true>>;

  void begin_end_update();

  iterator begin();

  iterator end();

  const_iterator cbegin();

  const_iterator cend();

  int64_t pointer_back_update();

  int64_t pointer_array_update();

  reverse_iterator rbegin();

  reverse_iterator rend();

  const_reverse_iterator crbegin();

  const_reverse_iterator crend();

  void insert(iterator iter, const T& value);

  template <typename... Args>
  void emplace(iterator iter, Args&&... args);

  void erase(iterator iter);

  T& operator[](size_t digit);

  T& at(size_t digit);

  const T& operator[](size_t digit) const;

  const T& at(size_t digit) const;

  bool empty() { return (size() == 0); }

  void delete_deq();

  Deque(int64_t count, const T& value, const Allocator& alloc = Allocator());

  Deque(int64_t count, const Allocator& alloc = Allocator());

  Deque(const Allocator& alloc = Allocator());

  Deque(const Deque& deque_old);

  Deque(Deque&& deque_old);

  Deque(std::initializer_list<T> init, const Allocator& alloc = Allocator());

  Deque<T, Allocator>& operator=(const Deque& deque_old);

  Deque<T, Allocator>& operator=(Deque&& deque_old);

  Allocator get_allocator() const { return alloc_; }

  void pointers_update();

  void push_back(const T& digit);

  void push_back(T&& digit);

  void push_front(T&& digit);

  void push_front(const T& digit);

  template <typename... Args>
  void emplace_back(Args&&... args);

  template <typename... Args>
  void emplace_front(Args&&... args);

  void push_back_fict();

  void swap(Deque& deque_old);

  void pop_back_fake();

  void pop_back();

  void pop_front_fake();

  void pop_front();

  T& back();

  T& front();

  size_t size() const;

  ~Deque();

 private:
  void relocate_push();

  Allocator alloc_;
  std::vector<T*> array_;
  int64_t pointer_back_array_ = 0;
  int64_t pointer_front_array_ = 1;
  int64_t pointer_back_ = kSize;
  int64_t pointer_front_ = -1;
};

template <typename T, typename Allocator>
void Deque<T, Allocator>::relocate_push() {
  try {
    std::vector<T*> new_array;
    if (array_.size() % 2 == 0) {
      new_array.resize(2 * array_.size());
    } else {
      new_array.resize(2 * array_.size() + 1);
    }
    for (size_t i = 0; i < new_array.size(); ++i) {
      T* temp = alloc_traits::allocate(alloc_, kSize);
      new_array[i] = temp;
    }
    for (size_t k = 0; k < array_.size(); ++k) {
      alloc_traits::deallocate(alloc_, new_array[((array_.size() + 1) / 2) + k],
                               kSize);
      new_array[((array_.size() + 1) / 2) + k] = array_[k];
    }
    pointer_back_array_ = pointer_back_array_ + ((array_.size() + 1) / 2);
    pointer_front_array_ = pointer_front_array_ + ((array_.size() + 1) / 2);
    array_ = new_array;
  } catch (...) {
    delete_deq();
  }
}

template <typename T, typename Allocator>
Deque<T, Allocator>::~Deque() {
  for (int64_t i = 0; i < static_cast<int64_t>(array_.size()); ++i) {
    if (pointer_back_array_ == pointer_front_array_) {
      for (int64_t j = pointer_back_; j <= pointer_front_; ++j) {
        alloc_traits::destroy(alloc_, array_[pointer_front_array_] + j);
      }
      break;
    }
    if (i == pointer_back_array_) {
      for (size_t j = pointer_back_; j < kSize; ++j) {
        alloc_traits::destroy(alloc_, array_[pointer_back_array_] + j);
      }
      continue;
    }
    if (i == pointer_front_array_) {
      for (int64_t j = 0; j <= pointer_front_; ++j) {
        alloc_traits::destroy(alloc_, array_[pointer_front_array_] + j);
      }
      continue;
    }
    if (i < pointer_front_array_ && i > pointer_back_array_) {
      for (size_t j = 0; j < kSize; ++j) {
        alloc_traits::destroy(alloc_, array_[i] + j);
      }
    }
  }
  for (size_t i = 0; i < array_.size(); ++i) {
    alloc_traits::deallocate(alloc_, array_[i], kSize);
  }
  array_.clear();
}

template <typename T, typename Allocator>
size_t Deque<T, Allocator>::size() const {
  if (pointer_back_ == kSize && pointer_front_ == -1) {
    return 0;
  }
  if (pointer_front_array_ == pointer_back_array_ && pointer_back_ != kSize &&
      pointer_front_ != -1) {
    return pointer_front_ - pointer_back_ + 1;
  }
  if (pointer_front_array_ == pointer_back_array_ && pointer_back_ == kSize) {
    return pointer_front_ + 1;
  }
  if (pointer_front_array_ == pointer_back_array_ && pointer_front_ != -1) {
    return kSize - pointer_back_;
  }
  return (pointer_front_array_ - pointer_back_array_ - 1) * kSize +
         (pointer_front_ - pointer_back_ + kSize) + 1;
}

template <typename T, typename Allocator>
T& Deque<T, Allocator>::front() {
  if (pointer_front_ == -1) {
    return array_[(array_.size() - 1) / 2][(kSize - 1)];
  }
  return array_[pointer_front_array_][pointer_front_];
}

template <typename T, typename Allocator>
T& Deque<T, Allocator>::back() {
  if (pointer_back_ == kSize) {
    return array_[(array_.size()) / 2][0];
  }
  return array_[pointer_back_array_][pointer_back_];
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::pop_front_fake() {
  if (pointer_front_ == pointer_back_ &&
      pointer_back_array_ == pointer_front_array_) {
    pointer_back_ = kSize;
    pointer_front_ = -1;
    return;
  }
  if (pointer_front_ == 0) {
    --pointer_front_array_;
    pointer_front_ = (kSize - 1);
    return;
  }
  --pointer_front_;
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::pop_front() {
  if (pointer_front_ == pointer_back_ &&
      pointer_back_array_ == pointer_front_array_) {
    alloc_traits::destroy(alloc_,
                          array_[pointer_front_array_] + pointer_front_);
    pointer_back_ = kSize;
    pointer_front_ = -1;
    return;
  }
  if (pointer_front_ == 0) {
    alloc_traits::destroy(alloc_,
                          array_[pointer_front_array_] + pointer_front_);
    --pointer_front_array_;
    pointer_front_ = (kSize - 1);
    return;
  }
  alloc_traits::destroy(alloc_, array_[pointer_front_array_] + pointer_front_);
  --pointer_front_;
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::pop_back_fake() {
  if (pointer_front_ == pointer_back_ &&
      pointer_back_array_ == pointer_front_array_) {
    pointer_back_ = kSize;
    pointer_front_ = -1;
    return;
  }
  if (pointer_back_ == (kSize - 1)) {
    ++pointer_back_array_;
    pointer_back_ = 0;
    return;
  }
  ++pointer_back_;
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::pop_back() {
  if (pointer_front_ == pointer_back_ &&
      pointer_back_array_ == pointer_front_array_) {
    alloc_traits::destroy(alloc_,
                          array_[pointer_front_array_] + pointer_front_);
    pointer_back_ = kSize;
    pointer_front_ = -1;
    return;
  }
  if (pointer_back_ == (kSize - 1)) {
    alloc_traits::destroy(alloc_, array_[pointer_back_array_] + pointer_back_);
    ++pointer_back_array_;
    pointer_back_ = 0;
    return;
  }
  alloc_traits::destroy(alloc_, array_[pointer_back_array_] + pointer_back_);
  ++pointer_back_;
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::push_front(const T& digit) {
  try {
    if (pointer_front_ == (kSize - 1) &&
        pointer_front_array_ == static_cast<int64_t>(array_.size()) - 1) {
      pointer_front_ = 0;
      relocate_push();
      ++pointer_front_array_;
      alloc_traits::construct(
          alloc_, array_[pointer_front_array_] + pointer_front_, digit);
      return;
    }
    if (pointer_front_ == (kSize - 1)) {
      pointer_front_ = 0;
      ++pointer_front_array_;
      alloc_traits::construct(
          alloc_, array_[pointer_front_array_] + pointer_front_, digit);
      return;
    }
    pointer_front_ = (++pointer_front_ % kSize);
    alloc_traits::construct(
        alloc_, array_[pointer_front_array_] + pointer_front_, digit);
    if (size() == 1) {
      pointers_update();
    }
  } catch (...) {
    pop_front_fake();
  }
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::push_front(T&& digit) {
  try {
    if (pointer_front_ == (kSize - 1) &&
        pointer_front_array_ == static_cast<int64_t>(array_.size()) - 1) {
      pointer_front_ = 0;
      relocate_push();
      ++pointer_front_array_;
      alloc_traits::construct(alloc_,
                              array_[pointer_front_array_] + pointer_front_,
                              std::move(digit));
      return;
    }
    if (pointer_front_ == (kSize - 1)) {
      pointer_front_ = 0;
      ++pointer_front_array_;
      alloc_traits::construct(alloc_,
                              array_[pointer_front_array_] + pointer_front_,
                              std::move(digit));
      return;
    }
    pointer_front_ = (++pointer_front_ % kSize);
    alloc_traits::construct(alloc_,
                            array_[pointer_front_array_] + pointer_front_,
                            std::move(digit));
    if (size() == 1) {
      pointers_update();
    }
  } catch (...) {
    pop_front_fake();
  }
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::push_back(const T& digit) {
  try {
    if (pointer_back_ == 0 && pointer_back_array_ == 0) {
      pointer_back_ = (kSize - 1);
      relocate_push();
      --pointer_back_array_;
      alloc_traits::construct(
          alloc_, array_[pointer_back_array_] + pointer_back_, digit);
      return;
    }
    if (pointer_back_ == 0) {
      pointer_back_ = (kSize - 1);
      --pointer_back_array_;
      alloc_traits::construct(
          alloc_, array_[pointer_back_array_] + pointer_back_, digit);
      return;
    }
    pointer_back_ = (--pointer_back_ % kSize);
    alloc_traits::construct(alloc_, array_[pointer_back_array_] + pointer_back_,
                            digit);
    if (size() == 1) {
      pointers_update();
    }
  } catch (...) {
    pop_back_fake();
    throw;
  }
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::push_back(T&& digit) {
  try {
    if (pointer_back_ == 0 && pointer_back_array_ == 0) {
      pointer_back_ = (kSize - 1);
      relocate_push();
      --pointer_back_array_;
      alloc_traits::construct(alloc_,
                              array_[pointer_back_array_] + pointer_back_,
                              std::move(digit));
      return;
    }
    if (pointer_back_ == 0) {
      pointer_back_ = (kSize - 1);
      --pointer_back_array_;
      alloc_traits::construct(alloc_,
                              array_[pointer_back_array_] + pointer_back_,
                              std::move(digit));
      return;
    }
    pointer_back_ = (--pointer_back_ % kSize);
    alloc_traits::construct(alloc_, array_[pointer_back_array_] + pointer_back_,
                            std::move(digit));
    if (size() == 1) {
      pointers_update();
    }
  } catch (...) {
    pop_back_fake();
    throw;
  }
}

template <typename T, typename Allocator>
template <typename... Args>
void Deque<T, Allocator>::emplace_front(Args&&... args) {
  try {
    if (pointer_front_ == (kSize - 1) &&
        pointer_front_array_ == static_cast<int64_t>(array_.size()) - 1) {
      pointer_front_ = 0;
      relocate_push();
      ++pointer_front_array_;
      alloc_traits::construct(alloc_,
                              array_[pointer_front_array_] + pointer_front_,
                              std::forward<Args>(args)...);
      return;
    }
    if (pointer_front_ == (kSize - 1)) {
      pointer_front_ = 0;
      ++pointer_front_array_;
      alloc_traits::construct(alloc_,
                              array_[pointer_front_array_] + pointer_front_,
                              std::forward<Args>(args)...);
      return;
    }
    pointer_front_ = (++pointer_front_ % kSize);
    alloc_traits::construct(alloc_,
                            array_[pointer_front_array_] + pointer_front_,
                            std::forward<Args>(args)...);
    if (size() == 1) {
      pointers_update();
    }
  } catch (...) {
    pop_front_fake();
  }
}

template <typename T, typename Allocator>
template <typename... Args>
void Deque<T, Allocator>::emplace_back(Args&&... args) {
  try {
    if (pointer_back_ == 0 && pointer_back_array_ == 0) {
      pointer_back_ = (kSize - 1);
      relocate_push();
      --pointer_back_array_;
      alloc_traits::construct(alloc_,
                              array_[pointer_back_array_] + pointer_back_,
                              std::forward<Args>(args)...);
      return;
    }
    if (pointer_back_ == 0) {
      pointer_back_ = (kSize - 1);
      --pointer_back_array_;
      alloc_traits::construct(alloc_,
                              array_[pointer_back_array_] + pointer_back_,
                              std::forward<Args>(args)...);
      return;
    }
    pointer_back_ = (--pointer_back_ % kSize);
    alloc_traits::construct(alloc_, array_[pointer_back_array_] + pointer_back_,
                            std::forward<Args>(args)...);
    if (size() == 1) {
      pointers_update();
    }
  } catch (...) {
    pop_back_fake();
    throw;
  }
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::push_back_fict() {
  try {
    if (pointer_back_ == 0 && pointer_back_array_ == 0) {
      pointer_back_ = (kSize - 1);
      relocate_push();
      --pointer_back_array_;
      alloc_traits::construct(alloc_,
                              array_[pointer_back_array_] + pointer_back_);
      return;
    }
    if (pointer_back_ == 0) {
      pointer_back_ = (kSize - 1);
      --pointer_back_array_;
      alloc_traits::construct(alloc_,
                              array_[pointer_back_array_] + pointer_back_);
      return;
    }
    pointer_back_ = (--pointer_back_ % kSize);
    alloc_traits::construct(alloc_,
                            array_[pointer_back_array_] + pointer_back_);
    if (size() == 1) {
      pointers_update();
    }
  } catch (...) {
    pop_back_fake();
    throw;
  }
}

template <typename T, typename Allocator>
Deque<T, Allocator>& Deque<T, Allocator>::operator=(const Deque& deque_old) {
  Deque<T, Allocator> copy = deque_old;
  if (alloc_traits::propagate_on_container_copy_assignment::value) {
    alloc_ = deque_old.alloc_;
  }
  swap(copy);
  return *this;
}

template <typename T, typename Allocator>
Deque<T, Allocator>& Deque<T, Allocator>::operator=(Deque&& deque_old) {
  Deque<T, Allocator> copy = std::move(deque_old);
  if (alloc_traits::propagate_on_container_copy_assignment::value) {
    alloc_ = deque_old.alloc_;
  }
  swap(copy);
  return *this;
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::swap(Deque& deque_old) {
  std::swap(pointer_back_, deque_old.pointer_back_);
  std::swap(pointer_front_, deque_old.pointer_front_);
  std::swap(pointer_back_array_, deque_old.pointer_back_array_);
  std::swap(pointer_front_array_, deque_old.pointer_front_array_);
  std::swap(array_, deque_old.array_);
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::pointers_update() {
  if (pointer_back_ == kSize) {
    pointer_back_ = pointer_front_;
    pointer_back_array_ = pointer_front_array_;
  }
  if (pointer_front_ == -1) {
    pointer_front_ = pointer_back_;
    pointer_front_array_ = pointer_back_array_;
  }
}

template <typename T, typename Allocator>
Deque<T, Allocator>::Deque(int64_t count, const T& value,
                           const Allocator& alloc) {
  alloc_ = alloc;
  for (int64_t i = 0; i < count / kSize + 1; i++) {
    T* temp = alloc_traits::allocate(alloc_, kSize);
    array_.push_back(temp);
  }
  try {
    for (int64_t i = 0; i < count; ++i) {
      push_back(value);
    }
  } catch (...) {
    delete_deq();
    throw;
  }
}

template <typename T, typename Allocator>
Deque<T, Allocator>::Deque(const Allocator& alloc) {
  try {
    alloc_ = alloc;
    for (int64_t i = 0; i < 2; i++) {
      T* temp = alloc_traits::allocate(alloc_, kSize);
      array_.push_back(temp);
    }
  } catch (...) {
    for (int64_t i = 0; i < 2; i++) {
      alloc_traits::deallocate(alloc_, array_[i], kSize);
    }
  }
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::delete_deq() {
  for (int64_t i = 0; i < static_cast<int64_t>(array_.size()); ++i) {
    if (pointer_back_array_ == pointer_front_array_) {
      for (int64_t j = pointer_back_; j <= pointer_front_; ++j) {
        alloc_traits::destroy(alloc_, array_[pointer_front_array_] + j);
      }
      break;
    }
    if (i == pointer_back_array_) {
      for (int64_t j = pointer_back_; j < kSize; ++j) {
        alloc_traits::destroy(alloc_, array_[pointer_back_array_] + j);
      }
      continue;
    }
    if (i == pointer_front_array_) {
      for (int64_t j = 0; j < pointer_front_; ++j) {
        alloc_traits::destroy(alloc_, array_[pointer_front_array_] + j);
      }
      continue;
    }
    if (i < pointer_front_array_ && i > pointer_back_array_) {
      for (int64_t j = 0; j < kSize; ++j) {
        alloc_traits::destroy(alloc_, array_[i] + j);
      }
    }
    alloc_traits::deallocate(alloc_, array_[i], kSize);
  }
  for (size_t i = 0; i < array_.size(); ++i) {
    alloc_traits::deallocate(alloc_, array_[i], kSize);
  }
  array_.clear();
}

template <typename T, typename Allocator>
Deque<T, Allocator>::Deque(int64_t count, const Allocator& alloc) {
  alloc_ = alloc;
  for (int64_t i = 0; i < count / kSize + 1; i++) {
    T* temp = alloc_traits::allocate(alloc_, kSize);
    array_.push_back(temp);
  }
  try {
    for (int64_t i = 0; i < count; ++i) {
      push_back_fict();
    }
  } catch (...) {
    delete_deq();
    throw;
  }
}

template <typename T, typename Allocator>
const T& Deque<T, Allocator>::at(size_t digit) const {
  int64_t new_digit = static_cast<int64_t>(digit);
  ++new_digit;
  if (digit >= size()) {
    throw std::out_of_range("ERROR");
  }
  if ((pointer_front_ - new_digit + 1) >= 0) {
    return array_[pointer_front_array_][pointer_front_ - new_digit + 1];
  }
  int64_t digit_help = new_digit;
  digit_help = digit_help - (pointer_front_ + 1);
  return array_[pointer_front_array_ - digit_help / kSize - 1]
               [(kSize - digit_help % kSize) % kSize];
}

template <typename T, typename Allocator>
const T& Deque<T, Allocator>::operator[](size_t digit) const {
  int64_t new_digit = static_cast<int64_t>(digit);
  ++new_digit;
  if ((pointer_front_ - new_digit + 1) >= 0) {
    return array_[pointer_front_array_][pointer_front_ - new_digit + 1];
  }
  int64_t digit_help = new_digit;
  digit_help = digit_help - (pointer_front_ + 1);
  return array_[pointer_front_array_ - digit_help / kSize - 1]
               [(kSize - digit_help % kSize) % kSize];
}

template <typename T, typename Allocator>
T& Deque<T, Allocator>::at(size_t digit) {
  int64_t new_digit = static_cast<int64_t>(digit);
  ++new_digit;
  if (digit >= size()) {
    throw std::out_of_range("ERROR");
  }
  if ((pointer_front_ - new_digit + 1) >= 0) {
    return array_[pointer_front_array_][pointer_front_ - new_digit + 1];
  }
  int64_t digit_help = new_digit;
  digit_help = digit_help - (pointer_front_ + 1);
  return array_[pointer_front_array_ - digit_help / kSize - 1]
               [(kSize - digit_help % kSize) % kSize];
}

template <typename T, typename Allocator>
T& Deque<T, Allocator>::operator[](size_t digit) {
  int64_t new_digit = static_cast<int64_t>(digit);
  ++new_digit;
  if ((pointer_front_ - new_digit + 1) >= 0) {
    return array_[pointer_front_array_][pointer_front_ - new_digit + 1];
  }
  int64_t digit_help = new_digit;
  digit_help = digit_help - (pointer_front_ + 1);
  return array_[pointer_front_array_ - digit_help / kSize - 1]
               [(kSize - digit_help % kSize) % kSize];
}

template <typename T, typename Allocator>
Deque<T, Allocator>::Deque(const Deque& deque_old) {
  try {
    alloc_ =
        alloc_traits::select_on_container_copy_construction(deque_old.alloc_);
    for (size_t i = 0; i < deque_old.array_.size(); i++) {
      T* temp = alloc_traits::allocate(alloc_, kSize);
      array_.push_back(temp);
    }
    for (int64_t i = 0; i < static_cast<int64_t>(array_.size()); ++i) {
      if (deque_old.pointer_back_array_ == deque_old.pointer_front_array_) {
        for (int64_t j = deque_old.pointer_back_; j <= deque_old.pointer_front_;
             ++j) {
          alloc_traits::construct(alloc_, array_[i] + j,
                                  deque_old.array_[i][j]);
          pointer_back_array_ = deque_old.pointer_back_array_;
          pointer_front_array_ = deque_old.pointer_front_array_;
          pointer_back_ = deque_old.pointer_back_;
          pointer_front_ = j;
        }
        break;
      }
      if (i == deque_old.pointer_back_array_) {
        for (int64_t j = deque_old.pointer_back_; j < kSize; ++j) {
          alloc_traits::construct(alloc_, array_[i] + j,
                                  deque_old.array_[i][j]);
        }
        continue;
      }
      if (i == deque_old.pointer_front_array_) {
        for (int64_t j = 0; j <= deque_old.pointer_front_; ++j) {
          alloc_traits::construct(alloc_, array_[i] + j,
                                  deque_old.array_[i][j]);
        }
        continue;
      }
      if (i < deque_old.pointer_front_array_ &&
          i > deque_old.pointer_back_array_) {
        for (int64_t j = 0; j < kSize; ++j) {
          alloc_traits::construct(alloc_, array_[i] + j,
                                  deque_old.array_[i][j]);
        }
      }
    }
    pointer_back_array_ = deque_old.pointer_back_array_;
    pointer_front_array_ = deque_old.pointer_front_array_;
    pointer_front_ = deque_old.pointer_front_;
    pointer_back_ = deque_old.pointer_back_;
  } catch (...) {
    delete_deq();
  }
}

template <typename T, typename Allocator>
Deque<T, Allocator>::Deque(std::initializer_list<T> init,
                           const Allocator& alloc) {
  try {
    alloc_ = alloc;
    for (size_t i = 0; i < init.size() / kSize + 1; ++i) {
      T* temp = alloc_traits::allocate(alloc_, kSize);
      array_.push_back(temp);
    }
    for (auto iter : init) {
      push_back(std::move(iter));
    }
  } catch (...) {
    delete_deq();
    throw;
  }
}

template <typename T, typename Allocator>
Deque<T, Allocator>::Deque(Deque&& deque_old) {
  try {
    pointer_front_ = deque_old.pointer_front_;
    pointer_back_ = deque_old.pointer_back_;
    pointer_back_array_ = deque_old.pointer_back_array_;
    pointer_front_array_ = deque_old.pointer_front_array_;
    alloc_ = std::move(deque_old.alloc_);
    array_ = std::move(deque_old.array_);
    deque_old.pointer_back_ = kSize;
    deque_old.pointer_front_ = -1;
    deque_old.pointer_front_array_ = 1;
    deque_old.pointer_back_array_ = 0;
  } catch (...) {
    delete_deq();
  }
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::erase(iterator iter) {
  int64_t place = iter.place();
  place--;
  for (int64_t i = place; i < static_cast<int64_t>(size()) - 1; ++i) {
    at(i) = at(i + 1);
  }
  pop_back();
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::insert(iterator iter, const T& value) {
  if (iter.set()) {
    push_back(value);
    return;
  }
  int64_t place = iter.place();
  place--;
  int64_t new_size = size();
  push_back(at(new_size - 1));
  for (int64_t i = new_size; i > place; --i) {
    at(i) = at(i - 1);
  }
  at(place) = value;
}

template <typename T, typename Allocator>
template <typename... Args>
void Deque<T, Allocator>::emplace(iterator iter, Args&&... args) {
  if (iter.set()) {
    emplace_back(std::forward<Args>(args)...);
    return;
  }
  int64_t place = iter.place();
  place--;
  int64_t new_size = size();
  push_back(at(new_size - 1));
  for (int64_t i = new_size; i > place; --i) {
    at(i) = at(i - 1);
  }
  at(place) = T(std::forward<Args>(args)...);
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::const_reverse_iterator
Deque<T, Allocator>::crend() {
  begin_end_update();
  return reverse_iterator(std::make_reverse_iterator(cbegin()));
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::const_reverse_iterator
Deque<T, Allocator>::crbegin() {
  begin_end_update();
  return reverse_iterator(std::make_reverse_iterator(cend()));
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::reverse_iterator Deque<T, Allocator>::rend() {
  begin_end_update();
  return reverse_iterator(std::make_reverse_iterator(begin()));
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::reverse_iterator Deque<T, Allocator>::rbegin() {
  begin_end_update();
  return reverse_iterator(std::make_reverse_iterator(end()));
}

template <typename T, typename Allocator>
int64_t Deque<T, Allocator>::pointer_array_update() {
  if (pointer_back_ == 0) {
    return (pointer_back_array_ - 1);
  }
  return (pointer_back_array_);
}

template <typename T, typename Allocator>
int64_t Deque<T, Allocator>::pointer_back_update() {
  if (pointer_back_ == 0) {
    return kSize - 1;
  }
  return (pointer_back_ - 1);
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::const_iterator Deque<T, Allocator>::cend() {
  begin_end_update();
  if (empty()) {
    return const_iterator(this);
  }
  return const_iterator(&array_[pointer_array_update()][pointer_back_update()],
                        this, pointer_back_update(), pointer_array_update());
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::const_iterator Deque<T, Allocator>::cbegin() {
  begin_end_update();
  if (empty()) {
    return const_iterator(this);
  }
  return const_iterator(&array_[pointer_front_array_][pointer_front_], this,
                        pointer_front_, pointer_front_array_);
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::iterator Deque<T, Allocator>::end() {
  begin_end_update();
  if (empty()) {
    return iterator(this);
  }
  return iterator(&array_[pointer_array_update()][pointer_back_update()], this,
                  pointer_back_update(), pointer_array_update());
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::iterator Deque<T, Allocator>::begin() {
  begin_end_update();
  if (empty()) {
    return iterator(this);
  }
  return iterator(&array_[pointer_front_array_][pointer_front_], this,
                  pointer_front_, pointer_front_array_);
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::begin_end_update() {
  if (pointer_back_ == 0 && pointer_back_array_ == 0) {
    relocate_push();
  }
  if (pointer_front_ == (kSize - 1) &&
      pointer_front_array_ == static_cast<int64_t>(array_.size() - 1)) {
    relocate_push();
  }
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template common_iterator<IsConst>::difference_type
Deque<T, Allocator>::common_iterator<IsConst>::operator-(
    const common_iterator& other) const {
  if (digit_pointer_ == 0 && other.digit_pointer_ == 0) {
    return 0;
  }
  if (digit_pointer_ == 0 && other.digit_pointer_ == 0) {
    return 0;
  }
  return (place() - other.place());
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template common_iterator<IsConst>
Deque<T, Allocator>::common_iterator<IsConst>::operator-(int64_t digit) const {
  common_iterator new_iter = *this;
  if (digit_pointer_ == 0 && digit == 0) {
    return new_iter;
  }
  if (digit_pointer_ == 0 && digit == 0) {
    return new_iter;
  }
  if (digit_pointer_ == -1 && pointer_array_ == -1) {
    return new_iter;
  }
  int64_t new_digit = (place() - digit);
  operators(new_digit, new_iter);
  new_iter.function_update();
  return new_iter;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template common_iterator<IsConst>
Deque<T, Allocator>::common_iterator<IsConst>::operator+(int64_t digit) const {
  common_iterator new_iter = *this;
  if (digit_pointer_ == 0 && digit == 0) {
    return new_iter;
  }
  if (digit_pointer_ == 0 && digit == 0) {
    return new_iter;
  }
  if (digit_pointer_ == -1 && pointer_array_ == -1) {
    return new_iter;
  }
  int64_t new_digit = place() + digit;
  operators(new_digit, new_iter);
  new_iter.function_update();
  return new_iter;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template common_iterator<IsConst>&
Deque<T, Allocator>::common_iterator<IsConst>::operator-=(int64_t digit) {
  int64_t new_digit = place() - digit;
  operators(new_digit, *this);
  function_update();
  return *this;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template common_iterator<IsConst>&
Deque<T, Allocator>::common_iterator<IsConst>::operator+=(int64_t digit) {
  int64_t new_digit = place() + digit;
  operators(new_digit, *this);
  function_update();
  return *this;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template common_iterator<IsConst>
Deque<T, Allocator>::common_iterator<IsConst>::operator++(int) {
  common_iterator new_iter = *this;
  int64_t new_digit = place() + 1;
  operators(new_digit, *this);
  function_update();
  return new_iter;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template common_iterator<IsConst>&
Deque<T, Allocator>::common_iterator<IsConst>::operator=(
    const common_iterator& other) {
  pointer_ = other.pointer_;
  deque_ = other.deque_;
  pointer_array_ = other.pointer_array_;
  digit_pointer_ = other.digit_pointer_;
  return *this;
}

template <typename T, typename Allocator>
template <bool IsConst>
void Deque<T, Allocator>::common_iterator<IsConst>::function_update() {
  if (digit_pointer_ == kSize) {
    digit_pointer_ = 0;
    ++pointer_array_;
    pointer_ = &((*deque_).array_[pointer_array_][digit_pointer_]);
  }
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template common_iterator<IsConst>
Deque<T, Allocator>::common_iterator<IsConst>::operator--(int) {
  common_iterator new_iter = *this;
  int64_t new_digit = place() - 1;
  operators(new_digit, *this);
  function_update();
  return new_iter;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template common_iterator<IsConst>&
Deque<T, Allocator>::common_iterator<IsConst>::operator--() {
  int64_t new_digit = place() - 1;
  operators(new_digit, *this);
  function_update();
  return *this;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template common_iterator<IsConst>&
Deque<T, Allocator>::common_iterator<IsConst>::operator++() {
  int64_t new_digit = place() + 1;
  operators(new_digit, *this);
  function_update();
  return *this;
}

template <typename T, typename Allocator>
template <bool IsConst>
void Deque<T, Allocator>::common_iterator<IsConst>::operators(
    int64_t& new_digit, common_iterator& iter) const {
  if ((iter.deque_)->pointer_front_ - new_digit + 1 >= 0) {
    iter.pointer_ =
        &(iter.deque_)
             ->array_[(iter.deque_)->pointer_front_array_]
                     [(iter.deque_)->pointer_front_ - new_digit + 1];
    iter.digit_pointer_ = (iter.deque_)->pointer_front_ - new_digit + 1;
    iter.pointer_array_ = (iter.deque_)->pointer_front_array_;
    return;
  }
  int64_t digit_help = static_cast<int64_t>(new_digit);
  digit_help = digit_help - ((iter.deque_)->pointer_front_ + 1);
  iter.pointer_ =
      &(iter.deque_)
           ->array_[(iter.deque_)->pointer_front_array_ - digit_help / kSize -
                    1][kSize - digit_help % kSize];
  iter.digit_pointer_ = kSize - digit_help % kSize;
  iter.pointer_array_ =
      (iter.deque_)->pointer_front_array_ - digit_help / kSize - 1;
}

template <typename T, typename Allocator>
template <bool IsConst>
void Deque<T, Allocator>::common_iterator<IsConst>::operators(
    int64_t& new_digit, common_iterator& iter) {
  if ((iter.deque_)->pointer_front_ - new_digit + 1 >= 0) {
    iter.pointer_ =
        &(iter.deque_)
             ->array_[(iter.deque_)->pointer_front_array_]
                     [(iter.deque_)->pointer_front_ - new_digit + 1];
    iter.digit_pointer_ = (iter.deque_)->pointer_front_ - new_digit + 1;
    iter.pointer_array_ = (iter.deque_)->pointer_front_array_;
    return;
  }
  int64_t digit_help = static_cast<int64_t>(new_digit);
  digit_help = digit_help - ((iter.deque_)->pointer_front_ + 1);
  iter.pointer_ =
      &(iter.deque_)
           ->array_[(iter.deque_)->pointer_front_array_ - digit_help / kSize -
                    1][kSize - digit_help % kSize];
  iter.digit_pointer_ = kSize - digit_help % kSize;
  iter.pointer_array_ =
      (iter.deque_)->pointer_front_array_ - digit_help / kSize - 1;
}

template <typename T, typename Allocator>
template <bool IsConst>
bool Deque<T, Allocator>::common_iterator<IsConst>::end_pointer() {
  if ((*deque_).pointer_back_ == (digit_pointer_ + 1) &&
      (*deque_).pointer_back_array_ == pointer_array_) {
    return true;
  }
  return (((*deque_).pointer_back_ == (kSize - 1) && digit_pointer_ == 0) &&
          ((*deque_).pointer_back_array_ == (pointer_array_ + 1)));
}

template <typename T, typename Allocator>
template <bool IsConst>
int64_t Deque<T, Allocator>::common_iterator<IsConst>::place() const {
  if (pointer_array_ == (*deque_).pointer_front_array_) {
    return (*deque_).pointer_front_ - digit_pointer_ + 1;
  }
  return ((*deque_).pointer_front_array_ - pointer_array_ - 1) * kSize +
         (kSize - digit_pointer_) + (*deque_).pointer_front_ + 1;
}

template <typename T, typename Allocator>
template <bool IsConst>
int64_t Deque<T, Allocator>::common_iterator<IsConst>::place() {
  if (pointer_array_ == (*deque_).pointer_front_array_) {
    return (*deque_).pointer_front_ - digit_pointer_ + 1;
  }
  return ((*deque_).pointer_front_array_ - pointer_array_ - 1) * kSize +
         (kSize - digit_pointer_) + (*deque_).pointer_front_ + 1;
}

template <typename T, typename Allocator>
template <bool IsConst>
Deque<T, Allocator>::common_iterator<IsConst>::common_iterator(
    T* pointer, Deque<T, Allocator>* deque, int64_t digit_pointer,
    int64_t pointer_array)
    : pointer_(pointer),
      deque_(deque),
      pointer_array_(pointer_array),
      digit_pointer_(digit_pointer) {
  if ((*deque_).empty()) {
    pointer_array_ = -1;
    digit_pointer_ = -1;
  }
}
