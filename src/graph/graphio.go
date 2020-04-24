package graph

import (
	"fmt"
	"io"
	"os"
	"strconv"
	"time"

	"github.com/biogo/hts/sam"
	"github.com/will-rowe/gfa"
	"github.com/will-rowe/groot/src/version"
)

// Store stores the GROOT graphs, using the graphID as the lookup key
type Store map[uint32]*GrootGraph

// SaveGraphAsGFA is a method to convert and save a GrootGraph in GFA format
func (GrootGraph *GrootGraph) SaveGraphAsGFA(fileName string, totalKmers int) (int, error) {
	// a flag to prevent dumping graphs which had no reads map
	graphUsed := false
	t := time.Now()
	stamp := fmt.Sprintf("variation graph created by groot (version %v) at: %v", version.GetVersion, t.Format("Mon Jan _2 15:04:05 2006"))
	msg := fmt.Sprintf("this graph is approximately weighted using k-mer frequencies from projected read sketches (total k-mers projected across all graphs: %d)", totalKmers)
	// create a GFA instance
	newGFA := gfa.NewGFA()
	_ = newGFA.AddVersion(1)
	newGFA.AddComment([]byte(stamp))
	newGFA.AddComment([]byte(msg))
	// transfer all the GrootGraphNode content to the GFA instance
	for _, node := range GrootGraph.SortedNodes {

		// some nodes will be marked for skipping after pruning, ignore these
		if node.Marked {
			continue
		}

		// record if this graph has had reads map
		if (graphUsed == false) && (node.KmerFreq > 0) {
			graphUsed = true
		}
		segID := strconv.FormatUint(node.SegmentID, 10)
		// create the segment
		seg, err := gfa.NewSegment([]byte(segID), []byte(node.Sequence))
		if err != nil {
			return 0, err
		}
		// the k-mer count corresponds to the node weight, which is its share of the k-mers from the projected sketches
		kmerCount := fmt.Sprintf("KC:i:%d", int((node.KmerFreq)))
		ofs, err := gfa.NewOptionalFields([]byte(kmerCount))
		if err != nil {
			return 0, err
		}
		seg.AddOptionalFields(ofs)
		seg.Add(newGFA)
		// create the links
		for _, outEdge := range node.OutEdges {
			toSeg := strconv.FormatUint(outEdge, 10)
			link, err := gfa.NewLink([]byte(segID), []byte("+"), []byte(toSeg), []byte("+"), []byte("0M"))
			if err != nil {
				return 0, err
			}
			link.Add(newGFA)
		}
	}
	// don't save the graph if no reads aligned
	if graphUsed == false {
		return 0, nil
	}
	// create the paths
	for pathID, pathName := range GrootGraph.Paths {
		// some paths won't have complete coverage, and have had their lengths set to 0 - ignore these paths
		if GrootGraph.Lengths[pathID] == 0 {
			continue
		}
		segments, overlaps := [][]byte{}, [][]byte{}
		for _, node := range GrootGraph.SortedNodes {

			// some nodes will be marked for skipping after pruning, ignore these
			if node.Marked {
				continue
			}
			for _, id := range node.PathIDs {
				if id == pathID {
					segment := strconv.FormatUint(node.SegmentID, 10) + "+"
					overlap := strconv.Itoa(len(node.Sequence)) + "M"
					segments = append(segments, []byte(segment))
					overlaps = append(overlaps, []byte(overlap))
					break
				}
			}
		}
		// add the path
		path, err := gfa.NewPath(pathName, segments, overlaps)
		if err != nil {
			return 0, err
		}
		path.Add(newGFA)
	}
	// create a gfaWriter and write the GFA instance
	outfile, err := os.Create(fileName)
	if err != nil {
		return 0, err
	}
	defer outfile.Close()
	writer, err := gfa.NewWriter(outfile, newGFA)
	if err != nil {
		return 0, err
	}
	err = newGFA.WriteGFAContent(writer)
	return 1, nil
}

// LoadGFA reads a GFA file into a GFA struct
func LoadGFA(fileName string) (*gfa.GFA, error) {
	// load the GFA file
	fh, err := os.Open(fileName)
	reader, err := gfa.NewReader(fh)
	if err != nil {
		return nil, fmt.Errorf("can't read gfa file: %v", err)
	}
	// collect the GFA instance
	myGFA := reader.CollectGFA()
	// read the file
	for {
		line, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("error reading line in gfa file: %v", err)
		}
		if err := line.Add(myGFA); err != nil {
			return nil, fmt.Errorf("error adding line to GFA instance: %v", err)
		}
	}
	return myGFA, nil
}

// GetSAMrefs is a method to convert all paths held in graphStore to sam.References
func (graphStore Store) GetSAMrefs() (map[int][]*sam.Reference, error) {
	references := make(map[int][]*sam.Reference)
	for graphID, grootGraph := range graphStore {
		references[int(graphID)] = make([]*sam.Reference, len(grootGraph.Paths))
		for pathID, path := range grootGraph.Paths {
			reference, err := sam.NewReference(string(path), "", "", grootGraph.Lengths[pathID], nil, nil)
			if err != nil {
				return nil, err
			}
			references[int(graphID)][pathID] = reference
		}
	}
	return references, nil
}
