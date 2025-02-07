#pragma once
#include <iomanip>
#include <iostream>
#include <streambuf>

// A custom buffer that discards everything
class NullBuffer : public std::streambuf {
protected:
    // Called when there's data to write into the buffer
    int overflow(int c) override {
        // Just pretend we handled the character successfully
        return c;
    }
};

// A custom ostream that uses NullBuffer
class NullStream : public std::ostream {
public:
    NullStream() : std::ostream(&m_sb) {}
private:
    NullBuffer m_sb;
};

// Depending on DEBUG_LOG, choose where debugStream points
#ifdef DEBUG_LOG
   // If DEBUG_LOG is defined, debugStream is std::cout
   inline std::ostream& debugStream = std::cout;
#else
   // If not, debugStream is a null stream that discards everything
   inline NullStream debugNullStream;          // We'll instantiate it once
   inline std::ostream& debugStream = debugNullStream; 
#endif
