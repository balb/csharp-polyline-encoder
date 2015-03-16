// C# translation of Mark Rambow's Java reimplementation of Mark McClures Javascript PolylineEncoder
// Mark Walters 2008 - mark[at]sol5.co.uk

using System.Collections.Generic;

public class Track {

    private List<Trackpoint> trackpoints = new List<Trackpoint>();

    public List<Trackpoint> getTrackpoints() {
        return this.trackpoints;
    }

    public void setTrackpoints(List<Trackpoint> trackpoints) {
        this.trackpoints = trackpoints;
    }

    public void addTrackpoint(Trackpoint trkpt) {
        this.trackpoints.Add(trkpt);
    }

}
