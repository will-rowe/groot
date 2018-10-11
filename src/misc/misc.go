// contains some misc helper functions etc. for GROOT
package misc

import (
	"encoding/binary"
	"encoding/gob"
	"errors"
	"log"
	"os"
	"strings"

	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

// a function to throw error to the log and exit the program
func ErrorCheck(msg error) {
	if msg != nil {
		log.Fatal("encountered error: ", msg)
	}
}

// a function to check for required flags
func CheckRequiredFlags(flags *pflag.FlagSet) error {
	requiredError := false
	flagName := ""

	flags.VisitAll(func(flag *pflag.Flag) {
		requiredAnnotation := flag.Annotations[cobra.BashCompOneRequiredFlag]
		if len(requiredAnnotation) == 0 {
			return
		}

		flagRequired := requiredAnnotation[0] == "true"

		if flagRequired && !flag.Changed {
			requiredError = true
			flagName = flag.Name
		}
	})

	if requiredError {
		return errors.New("Required flag `" + flagName + "` has not been set")
	}

	return nil
}

// StartLogging is a function to start the log...
func StartLogging(logFile string) *os.File {
	logPath := strings.Split(logFile, "/")
	joinedLogPath := strings.Join(logPath[:len(logPath)-1], "/")
	if len(logPath) > 1 {
		if _, err := os.Stat(joinedLogPath); os.IsNotExist(err) {
			if err := os.MkdirAll(joinedLogPath, 0700); err != nil {
				log.Fatal("can't create specified directory for log")
			}
		}
	}
	logFH, err := os.OpenFile(logFile, os.O_WRONLY|os.O_CREATE|os.O_APPEND, 0644)
	if err != nil {
		log.Fatal(err)
	}
	return logFH
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
	Version    string
	Ksize      int
	SigSize    int
	JSthresh   float64
	ReadLength int
	Containment bool
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
