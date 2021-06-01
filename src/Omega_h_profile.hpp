#ifndef OMEGA_H_STACK_HPP
#define OMEGA_H_STACK_HPP

#include <Omega_h_timer.hpp>
#include <Omega_h_filesystem.hpp>
#include <cstring>
#include <vector>
#include <memory>
#ifdef OMEGA_H_USE_KOKKOS
#include <Omega_h_kokkos.hpp>
#endif

namespace Omega_h {

class Comm;
typedef std::shared_ptr<Comm> CommPtr;

namespace profile {

struct Strings {
  std::vector<char> chars;
  std::size_t save(char const* str) {
    auto ret = chars.size();
    while (true) {
      chars.push_back(*str);
      if (*str == '\0') break;
      ++str;
    }
    return ret;
  }
  char const* get(std::size_t i) const { return &chars[i]; }
};

static constexpr std::size_t invalid = std::numeric_limits<std::size_t>::max();

struct Frame {
  std::size_t parent;
  std::size_t first_child;
  std::size_t last_child;
  std::size_t next_sibling;
  std::size_t name_ptr;
  Now start_time;
  double total_runtime;
  std::size_t number_of_calls;
};

struct History {
  std::vector<Frame> frames;
  std::size_t current_frame;
  std::size_t last_root;
  Strings names;
  Now start_time;
  bool do_percent;
  double chop;
  bool add_filename;
  CommPtr comm;
  History(CommPtr comm = nullptr, bool dopercent=false, double chop=0.0, bool add_filename=false);
  History(const History& h);
  inline const char* get_name(std::size_t frame) const {
    return names.get(frames[frame].name_ptr);
  }
  inline std::size_t find_child_of(
      std::size_t parent_index, char const* name) const {
    if (parent_index == invalid) return find_root(name);
    for (std::size_t child = frames[parent_index].first_child; child != invalid;
         child = frames[child].next_sibling) {
      if (0 == std::strcmp(get_name(child), name)) {
        return child;
      }
    }
    return invalid;
  }
  inline std::size_t create_child_of(
      std::size_t parent_index, char const* name) {
    if (parent_index == invalid) return create_root(name);
    auto index = frames.size();
    frames.push_back(Frame());
    auto& frame = frames.back();
    auto& parent_frame = frames[parent_index];
    frame.parent = parent_index;
    frame.first_child = invalid;
    frame.last_child = invalid;
    auto old_last = parent_frame.last_child;
    if (old_last != invalid) {
      frames[old_last].next_sibling = index;
    } else {
      parent_frame.first_child = index;
    }
    frame.next_sibling = invalid;
    parent_frame.last_child = index;
    frame.name_ptr = names.save(name);
    frame.total_runtime = 0.0;
    frame.number_of_calls = 0;
    return index;
  }
  inline std::size_t create_child_of_current(char const* name) {
    return create_child_of(current_frame, name);
  }
  inline std::size_t find_root(char const* name) const {
    if (frames.empty()) return invalid;
    for (std::size_t i = 0; i != invalid; i = frames[i].next_sibling) {
      if (frames[i].parent == invalid &&
          (0 == std::strcmp(get_name(i), name))) {
        return i;
      }
    }
    return invalid;
  }
  inline std::size_t create_root(char const* name) {
    auto index = frames.size();
    frames.push_back(Frame());
    auto& frame = frames.back();
    frame.parent = invalid;
    frame.first_child = invalid;
    frame.last_child = invalid;
    if (index != 0) {
      frames[last_root].next_sibling = index;
    }
    last_root = index;
    frame.next_sibling = invalid;
    frame.name_ptr = names.save(name);
    frame.total_runtime = 0.0;
    frame.number_of_calls = 0;
    return index;
  }
  inline std::size_t find(char const* name) {
    return find_child_of(current_frame, name);
  }
  inline std::size_t create(char const* name) {
    return create_child_of(current_frame, name);
  }
  inline std::size_t find_or_create(char const* name) {
    return find_or_create_child_of(current_frame, name);
  }
  inline std::size_t find_or_create_child_of(
      std::size_t parent_index, char const* name) {
    auto found = find_child_of(parent_index, name);
    if (found != invalid) return found;
    return create_child_of(parent_index, name);
  }
  inline std::size_t push(char const* name) {
    std::size_t id = find_or_create(name);
    current_frame = id;
    return id;
  }
  inline void pop() { current_frame = frames[current_frame].parent; }
  inline void start(char const* const name) {
    auto id = push(name);
    frames[id].number_of_calls += 1;
    frames[id].start_time = now();
  }
  inline double measure_runtime() { 
    auto current_time = now();
    auto current_runtime = current_time - frames[current_frame].start_time;
    return current_runtime;
  }
  inline double measure_total_runtime() { 
    return frames[current_frame].total_runtime + measure_runtime();
  }
  inline void stop() {
    frames[current_frame].total_runtime = measure_total_runtime();
    pop();
  }
  std::size_t first(std::size_t parent) const;
  std::size_t next(std::size_t sibling) const;
  std::size_t parent(std::size_t child) const;
  std::size_t pre_order_next(std::size_t frame) const;
  double time(std::size_t frame) const;
  std::size_t calls(std::size_t frame) const;
};

OMEGA_H_DLL extern History* global_singleton_history;

void simple_print(profile::History const& history);
History invert(History const& h);
void print_time_sorted(History const& h);
void print_top_down_and_bottom_up(History const& h, double total_runtime);
void print_top_sorted(History const& h, double total_runtime);

}  // namespace profile
}  // namespace Omega_h

namespace Omega_h {

inline void begin_code(char const* name, char const* file=0) {
#ifdef OMEGA_H_USE_KOKKOS
  Kokkos::Profiling::pushRegion(name);
#endif
  if (profile::global_singleton_history) {
    if (file == 0) {
      file = "Omega_h";
    }
    if (profile::global_singleton_history->add_filename) {
      std::string str = ::Omega_h::filesystem::path(file).filename().string()+"::"+std::string(name);
      profile::global_singleton_history->start(str.c_str());
    } else {
      profile::global_singleton_history->start(name);
    }
  }
}

inline double get_runtime () {
  double runtime = 0.0;
  if (profile::global_singleton_history) {
    runtime = profile::global_singleton_history->measure_total_runtime();
  }
  return runtime;
}

inline void end_code() {
#ifdef OMEGA_H_USE_KOKKOS
  Kokkos::Profiling::popRegion();
#endif
  if (profile::global_singleton_history) {
    profile::global_singleton_history->stop();
  }
}

struct ScopedTimer {
  ScopedTimer(char const* name, char const *file=0) { begin_code(name, file); }
  ~ScopedTimer() { end_code(); }
  ScopedTimer(ScopedTimer const&) = delete;
  ScopedTimer(ScopedTimer&&) = delete;
  ScopedTimer& operator=(ScopedTimer const&) = delete;
  ScopedTimer& operator=(ScopedTimer&&) = delete;
  inline double total_runtime() { return get_runtime(); }
};

}  // namespace Omega_h

#define OMEGA_H_TIME_FUNCTION                                                  \
  ::Omega_h::ScopedTimer omega_h_scoped_function_timer(__FUNCTION__, __FILE__)

#endif
