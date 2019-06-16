/**
@author Patrick Baus
@version 1.1.0 09/17/2018
*/
#ifndef COBS_CPP_H
#define COBS_CPP_H

#include <stdint.h>	 // uint8_t, etc.
#include <stddef.h>	 // size_t

namespace cobs {
   
    static const size_t encode(
        uint8_t* buffer, 
        const size_t size, 
        const size_t offset) __attribute__((unused));
    static const size_t encode(
        uint8_t* buffer, 
        const size_t size, 
        const size_t offset) {

      if (size > 254)
  	    return 0;

      uint8_t* startOfData = buffer + offset;
      uint8_t* endOfBlock = &startOfData[size-1];
      uint8_t* cursor;

      *(startOfData -1) = 0x00;

      do {
          for (cursor = endOfBlock; 
          *cursor != 0x00; cursor--) 
          {};

          *cursor = endOfBlock - cursor + 1;

          endOfBlock = cursor - 1;
      } while (cursor > startOfData);

      if (*(startOfData - 1) == 0x00) {
          *(startOfData - 1) = 0x01;
      }
      return size + 1;
    }

    static const size_t decode(
        uint8_t* buffer, 
        const size_t size) __attribute__((unused));
    static const size_t decode(
        uint8_t* buffer, 
        const size_t size) {
    
        if (size < 1 || size > 255)
            return 0;

        uint8_t tmp = 0;
        uint8_t* endOfBuffer = buffer + size;

        do {
            tmp = *buffer;   
            *buffer = 0x00;
            buffer += tmp;
        } while(buffer < endOfBuffer);

        return size -1;
    }
}
#endif
