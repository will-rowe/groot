// contains some misc helper functions etc. for GROOT
package misc

import (
	"encoding/binary"
	"encoding/gob"
	"log"
	"os"
)

// a function to throw error to the log and exit the program
func ErrorCheck(msg error) {
	if msg != nil {
		log.Fatal("encountered error: ", msg)
	}
}

// a function to print an array of unsigned integers as a string - taken from https://github.com/ekzhu/minhash-lsh TODO: benchmark with other stringify options
func Stringify(sig []uint64) string {
	hashValueSize := 8
	s := make([]byte, hashValueSize*len(sig))
	buf := make([]byte, 8)
	for i, v := range sig {
		binary.LittleEndian.PutUint64(buf, v)
		copy(s[i*hashValueSize:(i+1)*hashValueSize], buf[:hashValueSize])
	}
	return string(s)
}

/*
  A type to save the command information
*/
type IndexInfo struct {
	Ksize      int
	SigSize    int
	JSthresh   float64
	ReadLength int
}

// method to dump the info to file
func (self *IndexInfo) Dump(path string) error {
	file, err := os.Create(path)
	if err == nil {
		encoder := gob.NewEncoder(file)
		encoder.Encode(self)
	}
	file.Close()
	return err
}

// method to load info from file
func (self *IndexInfo) Load(path string) error {
	file, err := os.Open(path)
	if err == nil {
		decoder := gob.NewDecoder(file)
		err = decoder.Decode(self)
	}
	file.Close()
	return err
}
