#include <iostream>
#include <fstream>
#include <string>

class Logger {
private:
    std::ofstream logfile;
public:
    Logger() {}

    void open(const TString& filename) {
        // Close any previously opened file
        if (logfile.is_open()) {
            logfile.close();
        }
        logfile.open(filename);
        if (!logfile.is_open()) {
            std::cerr << "Error: Unable to open file: " << filename << std::endl;
        }
    }

  void close(){
    if (logfile.is_open()){
      logfile.close();
    }
  }

    template <typename T>
    Logger& operator<<(const T& data) {
      std::cout << data;    // Print to console
        if (logfile.is_open()) {
	  logfile << data;  // Write to file if file is open
 
        }
        return *this;
    }
   
};

