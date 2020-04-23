// contains some misc helper functions etc. for GROOT
package misc

import (
	"errors"
	"fmt"
	"log"
	"os"
	"runtime"
	"strings"

	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

// ErrorCheck is a function to throw error to the log and exit the program
func ErrorCheck(msg error) {
	if msg != nil {
		log.Fatalf("terminated\n\nERROR --> %v\n\n", msg)
	}
}

// CheckRequiredFlags is a function to check for required flags before running GROOT
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

// CheckSTDIN is a function to check that STDIN can be read
func CheckSTDIN() error {
	stat, err := os.Stdin.Stat()
	if err != nil {
		return fmt.Errorf("error with STDIN")
	}
	if (stat.Mode() & os.ModeNamedPipe) == 0 {
		return fmt.Errorf("no STDIN found")
	}
	return nil
}

// CheckDir is a function to check that a directory exists
func CheckDir(dir string) error {
	if dir == "" {
		return fmt.Errorf("no directory specified")
	}
	if _, err := os.Stat(dir); err != nil {
		if os.IsNotExist(err) {
			return fmt.Errorf("directory does not exist: %v", dir)
		}
		return fmt.Errorf("can't access adirectory (check permissions): %v", dir)
	}
	return nil
}

// CheckFile is a function to check that a file can be read
func CheckFile(file string) error {
	if _, err := os.Stat(file); err != nil {
		if os.IsNotExist(err) {
			return fmt.Errorf("file does not exist: %v", file)
		}
		return fmt.Errorf("can't access file (check permissions): %v", file)
	}
	return nil
}

// CheckExt is a function to check the extensions of a file
func CheckExt(file string, exts []string) error {
	splitFilename := strings.Split(file, ".")
	finalIdx := len(splitFilename) - 1
	if splitFilename[finalIdx] == "gz" {
		finalIdx--
	}
	err := fmt.Errorf("file does not have recognised extension: %v", file)
	for _, ext := range exts {
		if splitFilename[finalIdx] == ext {
			err = nil
			break
		}
	}
	return err
}

// Uint64SliceEqual returns true if two uint64[] are identical
func Uint64SliceEqual(a []uint64, b []uint64) bool {
	if len(a) != len(b) {
		return false
	}
	for i, v := range a {
		if v != b[i] {
			return false
		}
	}
	return true
}

// PrintMemUsage outputs the current, total and OS memory being used. As well as the number
// of garage collection cycles completed.
// lifted from: https://golangcode.com/print-the-current-memory-usage/
func PrintMemUsage() string {
	var m runtime.MemStats
	runtime.ReadMemStats(&m)
	// For info on each, see: https://golang.org/pkg/runtime/#MemStats
	return fmt.Sprintf("[ Heap Allocations: %vMb, OS Memory: %vMb, Num. GC cycles: %v ]", bToMb(m.HeapAlloc), bToMb(m.Sys), m.NumGC)
}

// bToMb converts bytes to megabytes
func bToMb(b uint64) uint64 {
	return b / 1024 / 1024
}
