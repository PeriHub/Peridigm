//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#ifndef PERIDIGM_LOGGING_HPP
#define PERIDIGM_LOGGING_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <functional>
#include <memory>
#include <mpi.h>
#include <Peridigm_Version.hpp>
#include <Sacado.hpp>

namespace PeridigmNS {

  enum class LogLevel {
      DEBUG,
      INFO,
      WARNING,
      ERROR
  };

  using namespace std;

  inline string executeCommand(const string& command) {
      array<char, 128> buffer;
      string result;
      shared_ptr<FILE> pipe(popen(command.c_str(), "r"), pclose);
      if (!pipe) {
          throw runtime_error("popen() failed!");
      }
      while (!feof(pipe.get())) {
          if (fgets(buffer.data(), 128, pipe.get()) != nullptr) {
              result += buffer.data();
          }
      }
      return result;
  }

  inline string getGitHash(const string& sourceDir) {
      string command = "cd " + sourceDir + " && git rev-parse HEAD";
      string gitHash = executeCommand(command);
      // Remove any trailing newline characters
      // gitHash.erase(std::remove(gitHash.begin(), gitHash.end(), '\n'), gitHash.end());
      return gitHash;
  }

  inline bool hasLocalChanges(const string& sourceDir) {
    string command = "cd " + sourceDir + " && git status --porcelain";
    string statusOutput = executeCommand(command);
    return !statusOutput.empty();
  }

  #ifdef CMAKE_SOURCE_DIR
    const string sourceDir = CMAKE_SOURCE_DIR;
    const string gitHash = getGitHash(sourceDir);
    const bool hasChanges = hasLocalChanges(sourceDir);
  #endif
  #ifdef CMAKE_BUILD_TYPE
    const string buildType = CMAKE_BUILD_TYPE;
  #endif
  // int mpi_size = 0;
  // bool mpiEnabled = false;
  // #ifdef HAVE_MPI
  //   MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  //   if (mpi_size>0){
  //     mpiEnabled = true;
  //   }
  // #endif

  inline void init_log_file() {
    time_t currentTime = time(nullptr);
    string datetime = asctime(localtime(&currentTime));
    ofstream logfile("log.txt", ios::out | ios::trunc);
    logfile << "\n-- Peridigm" << endl;
    logfile << "-- version " << PeridigmNS::Peridigm_Version() << endl;
    logfile << "-- git hash " << gitHash << endl;
    if (hasChanges) {
      logfile << "-- local changes detected" << endl;
    }
    logfile << "-- datetime " << datetime << endl;
    logfile.close();
  }

  using LogFunction = function<void(ofstream&, const string&)>;

  inline void log(const LogFunction& log_func, const string& message) {
    ofstream logfile("log.txt", ios_base::app);
    log_func(logfile, message);
  }

  inline void debug_log(ofstream& logfile, const string& message, const string& file, int line) {
    int mpi_id = 0;
    if(false){
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
      logfile << "DEBUG" << " - mpi_id " << mpi_id << " - " << file << ":" << line << " - " << message << endl;
    }
    else{
      logfile << "DEBUG" << " - " << file << ":" << line << " - " << message << endl;
    }
  }

  inline void info_log(ofstream& logfile, const string& message) {
    logfile << "INFO" << " - " << message << endl;
    cout << "INFO" << " - " << message << endl;
  }

  inline void warning_log(ofstream& logfile, const string& message, const string& file, int line) {
    int mpi_id = 0;
    if(false){
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
      logfile << "WARNING" << " - mpi_id " << mpi_id << " - " << file << ":" << line << " - " << message << endl;
    }
    else{
      logfile << "WARNING" << " - " << file << ":" << line << " - " << message << endl;
    }
  }

  inline void error_log(ofstream& logfile, const string& message, const string& file, int line) {
    int mpi_id = 0;
    if(false){
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
      logfile << "ERROR" << " - mpi_id " << mpi_id << " - " << file << ":" << line << " - " << message << endl;
    }
    else{
      logfile << "ERROR" << " - " << file << ":" << line << " - " << message << endl;
    }
  }

  inline void log(LogLevel level, const string& file, int line, const string& message) {
    switch (level) {
      case LogLevel::DEBUG:
        log(bind(&debug_log, placeholders::_1, message, file, line), message);
        break;
      case LogLevel::INFO:
        log(info_log, message);
        break;
      case LogLevel::WARNING:
        log(bind(&warning_log, placeholders::_1, message, file, line), message);
        break;
      case LogLevel::ERROR:
        log(bind(&error_log, placeholders::_1, message, file, line), message);
        break;
    }
  }

  inline void log(LogLevel level, const string& file, int line, const float message) {
    log(level, file, line, to_string(message));
  }

  inline void log(LogLevel level, const string& file, int line,  const double message) {
    log(level, file, line, to_string(message));
  }

  inline void log(LogLevel level, const string& file, int line,  const int message) {
    log(level, file, line, to_string(message));
  }

  template <typename ScalarT>
  string sacado_to_string(const ScalarT value) {
    ostringstream oss;
    oss << value;
    return oss.str();
  }

  template <typename ScalarT>
  inline void log(LogLevel level, const string& file, int line,  const ScalarT message) {
    log(level, file, line, sacado_to_string(message));
  }

  #define LOG(level, message) log(level, __FILE__, __LINE__, message)

  inline void testForTermination(bool test, const string& file, int line, const string& message) {

    ofstream logfile("log.txt", ios_base::app);

    if(test){
      if(false){
        int mpi_id = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
        logfile << "ERROR" << " - mpi_id " << mpi_id << " - " << file << ":" << line << " - " << message << endl;
        cout << "ERROR" << " - mpi_id " << mpi_id << " - " << file << ":" << line << " - " << message << endl;
      }
      else{
        logfile << "ERROR" << " - " << file << ":" << line << " - " << message << endl;
        cout << "ERROR" << " - " << file << ":" << line << " - " << message << endl;
      }
      terminate();
    }
  }

  #define TestForTermination(test, message) PeridigmNS::testForTermination(test, __FILE__, __LINE__, message)
  
}

#endif