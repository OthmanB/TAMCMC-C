#include <iostream>
#include <iomanip>
#include <libconfig.h++>

using namespace libconfig;

int main(){

  Config cfg;

  // Read the file. If there is an error, report it and exit.
  try
  {
    cfg.readFile("config.cfg");
  }
  catch(const FileIOException &fioex)
  {
    std::cerr << "I/O error while reading file." << std::endl;
    return(EXIT_FAILURE);
  }
  catch(const ParseException &pex)
  {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    return(EXIT_FAILURE);
  }
  
}