// C# translation of Mark Rambow's Java reimplementation of Mark McClures Javascript PolylineEncoder
// Mark Walters 2008 - mark[at]sol5.co.uk

using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Text;

public class PolylineEncoder {

	private int numLevels = 18;

	private int zoomFactor = 2;

	private double verySmall = 0.00001;

	private bool forceEndpoints = true;

	private double[] zoomLevelBreaks;

	private Dictionary<String, Double> bounds;

	// constructor
	public PolylineEncoder(int numLevels, int zoomFactor, double verySmall,
			bool forceEndpoints) {

		this.numLevels = numLevels;
		this.zoomFactor = zoomFactor;
		this.verySmall = verySmall;
		this.forceEndpoints = forceEndpoints;

		this.zoomLevelBreaks = new double[numLevels];

		for (int i = 0; i < numLevels; i++) {
			this.zoomLevelBreaks[i] = verySmall
					* Math.Pow(this.zoomFactor, numLevels - i - 1);
		}
	}

	public PolylineEncoder() {
		this.zoomLevelBreaks = new double[numLevels];

		for (int i = 0; i < numLevels; i++) {
			this.zoomLevelBreaks[i] = verySmall
					* Math.Pow(this.zoomFactor, numLevels - i - 1);
		}
	}

	public static void main(String[] args) {

		// initialize trackpoints
		// dirty hack just to show how it works :)

		Track trk = new Track();

		trk.addTrackpoint(new Trackpoint(52.29834, 8.94328));
		trk.addTrackpoint(new Trackpoint(52.29767, 8.93614));
		trk.addTrackpoint(new Trackpoint(52.29322, 8.93301));
		trk.addTrackpoint(new Trackpoint(52.28938, 8.93036));
		trk.addTrackpoint(new Trackpoint(52.27014, 8.97475));

		// encodeSignedNumber(floor1e5(coordinate)));
		Console.WriteLine(createEncodings(trk, 17, 1));
        

	}

	/**
	 * Douglas-Peucker algorithm, adapted for encoding
	 * 
	 * @return Hashtable [EncodedPoints;EncodedLevels]
	 * 
	 */
	public Dictionary<String, String> dpEncode(Track track) {
		int i, maxLoc = 0;
		Stack<int[]> stack = new Stack<int[]>();
		double[] dists = new double[track.getTrackpoints().Count];
		double maxDist, absMaxDist = 0.0, temp = 0.0;
		int[] current;
		String encodedPoints, encodedLevels;

		if (track.getTrackpoints().Count > 2) {
			int[] stackVal = new int[] { 0, (track.getTrackpoints().Count - 1) };
			stack.Push(stackVal);

			while (stack.Count > 0) {
				current = stack.Pop();
				maxDist = 0;

				for (i = current[0] + 1; i < current[1]; i++) {
					temp = this.distance(track.getTrackpoints()[i], track
							.getTrackpoints()[current[0]], track
							.getTrackpoints()[current[1]]);
					if (temp > maxDist) {
						maxDist = temp;
						maxLoc = i;
						if (maxDist > absMaxDist) {
							absMaxDist = maxDist;
						}
					}
				}
				if (maxDist > this.verySmall) {
					dists[maxLoc] = maxDist;
					int[] stackValCurMax = { current[0], maxLoc };
					stack.Push(stackValCurMax);
					int[] stackValMaxCur = { maxLoc, current[1] };
					stack.Push(stackValMaxCur);
				}
			}
		}

		// Console.WriteLine("createEncodings(" + track.getTrackpoints().Count
		// + "," + dists.length + ")");
		encodedPoints = createEncodings(track.getTrackpoints(), dists);
		// Console.WriteLine("encodedPoints \t\t: " + encodedPoints);
		// encodedPoints.replace("\\","\\\\");
        encodedPoints = encodedPoints.Replace("\\", "\\\\"); // replace(encodedPoints, "\\", "\\\\");
		//Console.WriteLine("encodedPoints slashy?\t\t: " + encodedPoints);

		encodedLevels = encodeLevels(track.getTrackpoints(), dists, absMaxDist);
		//Console.WriteLine("encodedLevels: " + encodedLevels);

		Dictionary<String, String> hm = new Dictionary<String, String>();
		hm.Add("encodedPoints", encodedPoints);
        hm.Add("encodedLevels", encodedLevels);
		return hm;

	}

    //public String replace(String s, String one, String another) {
    //    // In a string replace one substring with another
    //    if (s.Equals(""))
    //        return "";
    //    String res = "";
    //    int i = s.IndexOf(one, 0);
    //    int lastpos = 0;
    //    while (i != -1) {
    //        res += s.Substring(lastpos, i) + another;
    //        lastpos = i + one.Length();
    //        i = s.IndexOf(one, lastpos);
    //    }
    //    res += s.Substring(lastpos); // the rest
    //    return res;
    //}

	/**
	 * distance(p0, p1, p2) computes the distance between the point p0 and the
	 * segment [p1,p2]. This could probably be replaced with something that is a
	 * bit more numerically stable.
	 * 
	 * @param p0
	 * @param p1
	 * @param p2
	 * @return
	 */
	public double distance(Trackpoint p0, Trackpoint p1, Trackpoint p2) {
		double u, result = 0.0;

		if (p1.getLatDouble() == p2.getLatDouble()
				&& p1.getLonDouble() == p2.getLonDouble()) {
			result = Math.Sqrt(Math.Pow(p2.getLatDouble() - p0.getLatDouble(), 2)
					+ Math.Pow(p2.getLonDouble() - p0.getLonDouble(), 2));
		} else {
			u = ((p0.getLatDouble() - p1.getLatDouble())
					* (p2.getLatDouble() - p1.getLatDouble()) + (p0
					.getLonDouble() - p1.getLonDouble())
					* (p2.getLonDouble() - p1.getLonDouble()))
					/ (Math.Pow(p2.getLatDouble() - p1.getLatDouble(), 2) + Math
							.Pow(p2.getLonDouble() - p1.getLonDouble(), 2));

			if (u <= 0) {
				result = Math.Sqrt(Math.Pow(p0.getLatDouble() - p1.getLatDouble(),
						2)
						+ Math.Pow(p0.getLonDouble() - p1.getLonDouble(), 2));
			}
			if (u >= 1) {
				result = Math.Sqrt(Math.Pow(p0.getLatDouble() - p2.getLatDouble(),
						2)
						+ Math.Pow(p0.getLonDouble() - p2.getLonDouble(), 2));
			}
			if (0 < u && u < 1) {
				result = Math.Sqrt(Math.Pow(p0.getLatDouble() - p1.getLatDouble()
						- u * (p2.getLatDouble() - p1.getLatDouble()), 2)
						+ Math.Pow(p0.getLonDouble() - p1.getLonDouble() - u
								* (p2.getLonDouble() - p1.getLonDouble()), 2));
			}
		}
		return result;
	}

	/**
	 * @param points
	 *            set the points that should be encoded all points have to be in
	 *            the following form: Latitude, longitude\n
	 */
	public static Track pointsToTrack(String points) {
		Track trk = new Track();

        using (StringReader st = new StringReader(points))
        {
            string line = st.ReadLine();
            while (line != null)
            {
                String[] pointStrings = line.Split(", ".ToCharArray());
                trk.addTrackpoint(new Trackpoint(double.Parse(pointStrings[0]),
                    double.Parse(pointStrings[1])));
                line = st.ReadLine();
            }
        }

        return trk;
	}

    // kmlLineStringToTrack should split the points by " " rather than using a StringReader() to read by line

    ///**
    // * @param LineString
    // *            set the points that should be encoded all points have to be in
    // *            the following form: Longitude,Latitude,Altitude"_"...
    // */
    //public static Track kmlLineStringToTrack(String points) {
    //    Track trk = new Track();

    //    using (StringReader st = new StringReader(points))
    //    {
    //        string line = st.ReadLine();
    //        while (line != null)
    //        {
    //            String[] pointStrings = line.Split(",".ToCharArray());
    //            trk.addTrackpoint(new Trackpoint(double.Parse(pointStrings[1]),
    //                double.Parse(pointStrings[0]), double.Parse(pointStrings[2]))); 
    //            line = st.ReadLine();
    //        }
    //    }

    //    return trk;
    //}

	/**
	 * Goolge cant show Altitude, but its in some GPS/GPX Files
	 * Altitude will be ignored here so far
	 * @param points
	 * @return
	 */
	public static Track pointsAndAltitudeToTrack(String points) {
		Console.WriteLine("pointsAndAltitudeToTrack");
		Track trk = new Track();

        using (StringReader st = new StringReader(points))
        {
            string line = st.ReadLine();
            while (line != null)
            {
                String[] pointStrings = line.Split(",".ToCharArray());
                trk.addTrackpoint(new Trackpoint(double.Parse(pointStrings[1]),
                    double.Parse(pointStrings[0])));
                Console.WriteLine(double.Parse(pointStrings[1]).ToString() + ", "
                        + double.Parse(pointStrings[0]).ToString());
                line = st.ReadLine();
            }
        }

		return trk;
	}

	private static int floor1e5(double coordinate) {
		return (int) Math.Floor(coordinate * 1e5);
	}

	private static String encodeSignedNumber(int num) {
		int sgn_num = num << 1;
		if (num < 0) {
			sgn_num = ~(sgn_num);
		}
		return (encodeNumber(sgn_num));
	}

	private static String encodeNumber(int num) {

        StringBuilder encodeString = new StringBuilder();

		while (num >= 0x20) {
			int nextValue = (0x20 | (num & 0x1f)) + 63;
			encodeString.Append((char) (nextValue));
			num >>= 5;
		}

		num += 63;
		encodeString.Append((char) (num));

		return encodeString.ToString();
	}

	/**
	 * Now we can use the previous function to march down the list of points and
	 * encode the levels. Like createEncodings, we ignore points whose distance
	 * (in dists) is undefined.
	 */
	private String encodeLevels(List<Trackpoint> points, double[] dists,
			double absMaxDist) {
		int i;
        StringBuilder encoded_levels = new StringBuilder();

		if (this.forceEndpoints) {
			encoded_levels.Append(encodeNumber(this.numLevels - 1));
		} else {
			encoded_levels.Append(encodeNumber(this.numLevels
					- computeLevel(absMaxDist) - 1));
		}
		for (i = 1; i < points.Count - 1; i++) {
			if (dists[i] != 0) {
				encoded_levels.Append(encodeNumber(this.numLevels
						- computeLevel(dists[i]) - 1));
			}
		}
		if (this.forceEndpoints) {
			encoded_levels.Append(encodeNumber(this.numLevels - 1));
		} else {
			encoded_levels.Append(encodeNumber(this.numLevels
					- computeLevel(absMaxDist) - 1));
		}
//		Console.WriteLine("encodedLevels: " + encoded_levels);
		return encoded_levels.ToString();
	}

	/**
	 * This computes the appropriate zoom level of a point in terms of it's
	 * distance from the relevant segment in the DP algorithm. Could be done in
	 * terms of a logarithm, but this approach makes it a bit easier to ensure
	 * that the level is not too large.
	 */
	private int computeLevel(double absMaxDist) {
		int lev = 0;
		if (absMaxDist > this.verySmall) {
			lev = 0;
			while (absMaxDist < this.zoomLevelBreaks[lev]) {
				lev++;
			}
			return lev;
		}
		return lev;
	}

	private String createEncodings(List<Trackpoint> points, double[] dists) {
        StringBuilder encodedPoints = new StringBuilder();

		double maxlat = 0, minlat = 0, maxlon = 0, minlon = 0;

		int plat = 0;
		int plng = 0;
		
		for (int i = 0; i < points.Count; i++) {

			// determin bounds (max/min lat/lon)
			if (i == 0) {
				maxlat = minlat = points[i].getLatDouble();
				maxlon = minlon = points[i].getLonDouble();
			} else {
				if (points[i].getLatDouble() > maxlat) {
					maxlat = points[i].getLatDouble();
				} else if (points[i].getLatDouble() < minlat) {
					minlat = points[i].getLatDouble();
				} else if (points[i].getLonDouble() > maxlon) {
					maxlon = points[i].getLonDouble();
				} else if (points[i].getLonDouble() < minlon) {
					minlon = points[i].getLonDouble();
				}
			}

			if (dists[i] != 0 || i == 0 || i == points.Count - 1) {
				Trackpoint point = points[i];

				int late5 = floor1e5(point.getLatDouble());
				int lnge5 = floor1e5(point.getLonDouble());

				int dlat = late5 - plat;
				int dlng = lnge5 - plng;

				plat = late5;
				plng = lnge5;

				encodedPoints.Append(encodeSignedNumber(dlat));
				encodedPoints.Append(encodeSignedNumber(dlng));

			}
		}

		Dictionary<String, Double> bounds = new Dictionary<String, Double>();
		bounds.Add("maxlat", maxlat);
        bounds.Add("minlat", minlat);
        bounds.Add("maxlon", maxlon);
        bounds.Add("minlon", minlon);

		this.setBounds(bounds);
		return encodedPoints.ToString();
	}

	private void setBounds(Dictionary<String, Double> bounds) {
		this.bounds = bounds;
	}

	public static Dictionary<String, String> createEncodings(Track track, int level, int step) {

		Dictionary<String, String> resultMap = new Dictionary<String, String>();
		StringBuilder encodedPoints = new StringBuilder();
		StringBuilder encodedLevels = new StringBuilder();

		List<Trackpoint> trackpointList = track.getTrackpoints();

		int plat = 0;
		int plng = 0;
		int counter = 0;

		int listSize = trackpointList.Count;

		Trackpoint trackpoint;

		for (int i = 0; i < listSize; i += step) {
			counter++;
			trackpoint = (Trackpoint) trackpointList[i];

			int late5 = floor1e5(trackpoint.getLatDouble());
			int lnge5 = floor1e5(trackpoint.getLonDouble());

			int dlat = late5 - plat;
			int dlng = lnge5 - plng;

			plat = late5;
			plng = lnge5;

			encodedPoints.Append(encodeSignedNumber(dlat)).Append(
					encodeSignedNumber(dlng));
			encodedLevels.Append(encodeNumber(level));

		}

		Console.WriteLine("listSize: " + listSize + " step: " + step
				+ " counter: " + counter);

		resultMap.Add("encodedPoints", encodedPoints.ToString());
		resultMap.Add("encodedLevels", encodedLevels.ToString());

		return resultMap;
	}

	public Dictionary<String, Double> getBounds() {
		return bounds;
	}
}
