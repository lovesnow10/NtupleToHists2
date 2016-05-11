#include "HplusRun.hpp"

int main(int argc, char **argv) {
  HplusRun *app = new HplusRun(argc, argv);
  app->initialize();
  app->run();
  app->finalize();
  return 0;
}
