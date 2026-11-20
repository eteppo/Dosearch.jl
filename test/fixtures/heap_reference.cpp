// Reference pop order of std::priority_queue (GCC libstdc++) with tied scores.
// Regenerate heap_ops.txt with: g++ -O1 -o heap heap_reference.cpp && ./heap > heap_ops.txt
#include <cstdio>
#include <queue>
#include <vector>
struct item { int score, id; };
struct comp { bool operator()(item const* a, item const* b) { return a->score < b->score; } };
int main() {
  unsigned s = 12345;
  auto rnd = [&]() { s = s * 1103515245u + 12345u; return (s >> 16) & 0x7fff; };
  std::vector<item> items(400);
  std::priority_queue<item*, std::vector<item*>, comp> q;
  int next = 0;
  for (int step = 0; step < 600; step++) {
    if (next < 400 && (q.empty() || rnd() % 3 != 0)) {
      items[next] = {int(rnd() % 5) - 2, next + 1};
      q.push(&items[next]);
      printf("push %d %d\n", items[next].id, items[next].score);
      next++;
    } else if (!q.empty()) {
      printf("pop %d\n", q.top()->id);
      q.pop();
    }
  }
  return 0;
}
