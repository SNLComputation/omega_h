#include <Omega_h_profile.hpp>
#include <Omega_h_print.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_dbg.hpp>
#include <algorithm>
#include <iostream>
#include <algorithm>
#include <queue>
#include <limits>
#include <iomanip>
#include <map>
#include <utility>

namespace Omega_h {
namespace profile {

OMEGA_H_DLL History* global_singleton_history = nullptr;

History::History(CommPtr comm_in, bool dopercent, double chop_in, bool add_filename_in) : 
  current_frame(invalid), last_root(invalid), start_time(now()), 
  do_percent(dopercent), chop(chop_in), add_filename(add_filename_in), comm(comm_in) {}

History::History(const History& h) {
  start_time = h.start_time;
  do_percent = h.do_percent;
  chop = h.chop;
  add_filename = h.add_filename;
  comm = h.comm;
}

std::size_t History::first(std::size_t parent_index) const {
  if (parent_index != invalid) return frames[parent_index].first_child;
  if (!frames.empty()) return 0;
  return invalid;
}

std::size_t History::next(std::size_t sibling) const {
  return frames[sibling].next_sibling;
}

std::size_t History::parent(std::size_t child) const {
  return frames[child].parent;
}

std::size_t History::pre_order_next(std::size_t frame) const {
  auto first2 = first(frame);
  if (first2 != invalid) {
    return first2;
  }
  auto next2 = next(frame);
  if (next2 != invalid) {
    return next2;
  }
  for (; frame != invalid; frame = parent(frame)) {
    auto next3 = next(frame);
    if (next3 != invalid) {
      return next3;
    }
  }
  return invalid;
}

double History::time(std::size_t frame) const {
  return frames[frame].total_runtime;
}

std::size_t History::calls(std::size_t frame) const {
  return frames[frame].number_of_calls;
}

struct PreOrderIterator {
  using reference = std::size_t;
  PreOrderIterator& operator++() {
    frame = history.pre_order_next(frame);
    return *this;
  }
  reference operator*() { return frame; }
  PreOrderIterator(History const& history_in, std::size_t i_in)
      : history(history_in), frame(i_in) {}
  bool operator!=(PreOrderIterator const& other) {
    return frame != other.frame;
  }
  History const& history;
  std::size_t frame;
};

static PreOrderIterator begin(History const& history) {
  return PreOrderIterator(history, history.frames.empty() ? invalid : 0);
}

static PreOrderIterator end(History const& history) {
  return PreOrderIterator(history, invalid);
}

static std::vector<std::size_t> compute_depths(History const& history) {
  std::vector<std::size_t> depths(history.frames.size());
  for (auto frame : history) {
    auto parent = history.parent(frame);
    if (parent != invalid) {
      depths[frame] = depths[parent] + 1;
    } else {
      depths[frame] = 0;
    }
  }
  return depths;
}

static void simple_print(
    History const& history, std::vector<std::size_t> const& depths) {
  for (auto frame : history) {
    std::size_t depth = depths[frame];
    for (std::size_t i = 0; i < depth; ++i) std::cout << "  ";
    std::cout << history.get_name(frame) << ' '
              << history.frames[frame].total_runtime << '\n';
  }
}

void simple_print(History const& history) {
  auto depths = compute_depths(history);
  simple_print(history, depths);
}

History invert(History const& h) {
  History invh(h);
  std::queue<std::size_t> q;
  for (std::size_t s = h.first(invalid); s != invalid; s = h.next(s)) {
    q.push(s);
  }
  while (!q.empty()) {
    auto node = q.front();
    q.pop();
    auto self_time = h.time(node);
    auto calls = h.calls(node);
    for (auto child = h.first(node); child != invalid; child = h.next(child)) {
      self_time -= h.time(child);
      q.push(child);
    }
    self_time = std::max(self_time,
        0.);  // floating-point may give negative epsilon instead of zero
    auto inv_node = invalid;
    for (; node != invalid; node = h.parent(node)) {
      auto name = h.get_name(node);
      inv_node = invh.find_or_create_child_of(inv_node, name);
      invh.frames[inv_node].total_runtime += self_time;
      invh.frames[inv_node].number_of_calls += calls;
    }
  }
  return invh;
}

static void print_time_sorted_recursive(History const& h, std::size_t frame,
    std::vector<std::size_t> const& depths, double total_runtime) {
  std::string percent = " ";
  double scale = 1.0;
  if (h.do_percent) {
    percent = "% ";
    scale = 100.0/total_runtime;
  }
  std::vector<std::size_t> child_frames;
  for (std::size_t child = h.first(frame); child != invalid;
       child = h.next(child)) {
    child_frames.push_back(child);
  }
  std::stable_sort(begin(child_frames), end(child_frames),
      [&](std::size_t a, std::size_t b) { return h.time(a) > h.time(b); });
  for (auto child : child_frames) {
    std::size_t depth = depths[child];
    if (h.time(child)*100.0/total_runtime >= h.chop) {
      for (std::size_t i = 0; i < depth; ++i) std::cout << "|  ";
      std::cout << h.get_name(child) << ' ' << h.time(child)*scale << percent 
                << h.calls(child) << '\n';
    }
    print_time_sorted_recursive(h, child, depths, total_runtime);
  }
}

void print_time_sorted(History const& h, double total_runtime) {
  auto depths = compute_depths(h);
  print_time_sorted_recursive(h, invalid, depths, total_runtime);
}

enum { TOP_AVE,
       TOP_MIN,
       TOP_MAX,
       TOP_NUM_ENTRIES };

static void gather(History const& h, std::size_t frame, std::map<std::string, std::vector<double>>& result) {
  for (std::size_t child = h.first(frame); child != invalid;
       child = h.next(child)) {
    if (result[h.get_name(child)].size() == 0) {
      result[h.get_name(child)].resize(TOP_NUM_ENTRIES);
    }
    result[h.get_name(child)][TOP_AVE] = h.time(child);
    result[h.get_name(child)][TOP_MIN] = h.time(child);
    result[h.get_name(child)][TOP_MAX] = h.time(child);
  }
}

void print_top_down_and_bottom_up(History const& h, double total_runtime) {
  auto coutflags( std::cout.flags() );
  if (h.do_percent) {
    std::cout << std::setprecision(2) << std::fixed;
  }
  std::cout << "\n";
  std::cout << "TOP-DOWN:\n";
  std::cout << "=========\n";
  print_time_sorted(h, total_runtime);
  auto h_inv = invert(h);
  std::cout << "\n";
  std::cout << "BOTTOM-UP:\n";
  std::cout << "==========\n";
  print_time_sorted(h_inv, total_runtime);
  std::cout.flags(coutflags);
}

void split_char_vec(const std::vector<char>& cvec,  std::vector<std::string>& res) {
  res.clear();
  size_t j = 0;
  const char* p = cvec.data();
  while(j < cvec.size()) {
    std::string s(p);
    res.push_back(s);
    int t = strlen(p);
    p += t + 1;
    j += t + 1;
  }
}

void sendrecv(History const& h, std::map<std::string, std::vector<double> >& result) {
  if (h.comm.get()) {
    std::vector<char> cvec;
    std::vector<double> dvec;
    if (h.comm->rank()) {
      for (auto i : result) {
        cvec.insert(cvec.end(), i.first.c_str(), i.first.c_str()+i.first.length()+1);
        dvec.push_back(i.second[TOP_AVE]);
      }
      h.comm->send(0, cvec);
      h.comm->send(0, dvec);
    } else {
      for (int irank = 1; irank < h.comm->size(); ++irank) {
        cvec.clear();
        dvec.clear();
        h.comm->recv(irank, cvec);
        h.comm->recv(irank, dvec);
        std::vector<std::string> res;
        split_char_vec(cvec, res);
        OMEGA_H_CHECK_OP(res.size(), ==, dvec.size());
        for (size_t i = 0; i < res.size(); ++i) {
          result[res[i]][TOP_AVE] += dvec[i];
          result[res[i]][TOP_MIN] = std::min(result[res[i]][TOP_MIN], dvec[i]);
          result[res[i]][TOP_MAX] = std::max(result[res[i]][TOP_MAX], dvec[i]);
        }
      }
    }
  }
}

void print_top_sorted(History const& h_in, double total_runtime) {
  auto h = invert(h_in);
  auto coutflags( std::cout.flags() );
  if (h.do_percent) {
    std::cout << std::setprecision(2) << std::fixed;
  }
  double sz =  h.comm->size();
  double total_runtime_ave = h.comm->allreduce(total_runtime, OMEGA_H_SUM) / sz;
  double total_runtime_min = h.comm->allreduce(total_runtime, OMEGA_H_MIN);
  double total_runtime_max = h.comm->allreduce(total_runtime, OMEGA_H_MAX);
  TASK_0_cout << "\n";
  TASK_0_cout << "TOP FUNCTIONS (self time, average of all ranks):\n";
  TASK_0_cout << "=============\n";
  std::map<std::string, std::vector<double>> result;
  gather(h, invalid, result);
  sendrecv(h, result);
  typedef std::pair<std::string, std::vector<double>> my_pair;
  std::vector<my_pair> sorted_result;
  double sum = 0.0;
  for (auto& i : result) {
    i.second[TOP_AVE] /= sz;
    sum += i.second[TOP_AVE];
    sorted_result.push_back(std::make_pair(i.first, i.second));
  }
  TASK_0_cout << "total_runtime= " << total_runtime_ave << " [s] monitored functions= " << sum 
              << " [s] unmonitored= " << 100.0*(total_runtime_ave - sum)/total_runtime_ave << "%" << std::endl;
  std::vector<double> vv = {(total_runtime_ave-sum), (total_runtime_min-sum), (total_runtime_max-sum)};
  sorted_result.push_back(std::make_pair("unmonitored functions", vv));
  std::stable_sort(sorted_result.begin(), sorted_result.end(),
            [](const my_pair& a, const my_pair& b) -> bool
            { 
              return a.second[0] > b.second[0];
            });  

  auto cflags( std::cout.flags() );
  std::string percent = "  ";
  std::string ul = "  ";
  double scale = 1.0;
  int width = 14;
  if (h.do_percent) {
    percent = "% ";
    ul = "- ";
    scale = 100.0/total_runtime_ave;
    width = 8;
  }
  TASK_0_cout << std::right 
              << std::setw(width) << "Ave" << percent 
              << std::setw(width) << "Min" << percent
              << std::setw(width) << "Max" << percent
              << std::setw(width) << "Max/Min" << "  "
              << " Name"
              << std::endl;
  TASK_0_cout << std::right 
              << std::setw(width) << "---" << ul
              << std::setw(width) << "---" << ul
              << std::setw(width) << "---" << ul
              << std::setw(width) << "-------" << "  "
              << " ----"
              << std::endl;
  std::cout.flags(cflags);

  for (auto i : sorted_result) {
    auto cflags( std::cout.flags() );
    double val = i.second[TOP_AVE];
    if (val*100.0/total_runtime_ave >= h.chop) {
      TASK_0_cout << std::right
                  << std::setw(width) << val*scale << percent 
                  << std::setw(width) << i.second[TOP_MIN]*scale << percent
                  << std::setw(width) << i.second[TOP_MAX]*scale << percent
                  << std::setw(width) << i.second[TOP_MAX]/i.second[TOP_MIN] << "   ";
      std::cout.flags(cflags);
      TASK_0_cout << i.first << std::endl;
    }
    std::cout.flags(coutflags);
  }
}

}  // namespace profile
}  // namespace Omega_h
