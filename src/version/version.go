package version

import "fmt"

// Version components are set at build time via ldflags.
var (
	major = "1"
	minor = "1"
	patch = "2"
)

// GetVersion returns the full version string for the current GROOT software
func GetVersion() string {
	return fmt.Sprintf("%s.%s.%s", major, minor, patch)
}

// GetBaseVersion returns the major minor version string for the current GROOT software
func GetBaseVersion() string {
	return fmt.Sprintf("%s.%s", major, minor)
}
