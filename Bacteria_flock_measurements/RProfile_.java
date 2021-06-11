/**
 * Plugin for extraction of exudate activity / quantity along a root
 * Lionel Dupuy 2016
 * JHI - Plant Systems Modelling
 */
 
import ij.*;
import ij.measure.*;
import ij.plugin.*;
import ij.plugin.filter.*;
import ij.process.*;
import ij.process.ImageProcessor.*;
import ij.process.AutoThresholder.Method;
import ij.gui.*;
import ij.gui.PolygonRoi;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.ContrastEnhancer;

import java.awt.*;
//import java.lang.Math.*;
import ij.plugin.frame.RoiManager;
import ij.plugin.filter.MaximumFinder;
import ij.process.ImageConverter;
import ij.process.EllipseFitter;
import ij.plugin.ImageCalculator;

import ij.gui.Roi;
import ij.gui.Roi.*;
import ij.gui.ShapeRoi;
import ij.plugin.frame.RoiManager;
import java.util.*;
import java.awt.font.*;

import java.io.*;

public class RProfile_ implements PlugIn {
	@SuppressWarnings("unchecked") 
	/********************************************************************************************
	 * Parameters of the profile algorithm
	 * Change here for best fit to your image
	 */
	int w_profile = 200; 				// fwidth of the root profile in pixels
	int l_profile = 10;//(int)(w_profile/2);		// length of root segments in pixels (after resampling)	
	/**
	 *********************************************************************************************
	 */
	
	// constants
	int n_slice;
	ImagePlus imp;
	// the tracing
	float[] x; // x coordinates of the tracing 
	float[] y; // y coordinate of the tracing
	float[] L; // distance from the tip (curvilinear abscisa)
	float[] NX; // normal vector coord x
	float[] NY; // normal vector coord y
	
	// float image
	float[][] I;

	
	ImageStack main_stack;
	int width = 0;
	int height = 0;
	public void run(String arg) 
	{
		if (IJ.versionLessThan("1.26i"))
			return;
		//IJ.run("32-bit");
		imp = IJ.getImage();
		
		// Image Data
		ImageProcessor IP = imp.getProcessor();
		width = IP.getWidth();
		height = IP.getHeight();
		I = new float[width][height];
		I = IP.getFloatArray();

		double[][] profile =  Read_Profile_From_ROI();
		profile =  Refine_Profile(profile);
		Export_Data(profile);
		Refresh_Roi();
	}


	/******************************************************************************************
	 * Read ROI data
	 */	

	private double[][] Read_Profile_From_ROI()
	{
		double[][] profile;
		Roi roi = imp.getRoi();
		PolygonRoi line_roi = (PolygonRoi)roi;
		FloatPolygon Poly = line_roi.getFloatPolygon();
		x = Poly.xpoints; 
		y = Poly.ypoints; 
		L = new float[x.length];
		NX = new float[x.length];
		NY = new float[x.length];		
		int n = line_roi.getNCoordinates(); 
		if (n>2) { 
			profile = new double[n][w_profile];
			for (int i=0;i<n;i++)
			{
				// determine normal vector
				float[] N = {0,0};
				float[] T = {0,0};
				
				if (i==0)
				{
					T = get_tangent(x[i],y[i],x[i+1],y[i+1]);
					N[0] = -T[1]; N[1] = T[0];
					// update
					L[i] = 0;
					NX[i] = N[0];
					NY[i] = N[1];
				}
				else if (i==n-1)
				{
					T = get_tangent(x[i-1],y[i-1],x[i],y[i]);
					N[0] = -T[1]; N[1] = T[0];	
					// update
					L[i] = L[i-1] + (float)Math.sqrt(Math.pow(x[i]-x[i-1],2) + Math.pow(y[i]-y[i-1],2));
					NX[i] = N[0];
					NY[i] = N[1];
				}
				else
				{
					T = get_tangent(x[i-1],y[i-1],x[i+1],y[i+1]);
					N[0] = -T[1]; N[1] = T[0];
					// update
					L[i] = L[i-1] + (float)Math.sqrt(Math.pow(x[i]-x[i-1],2) + Math.pow(y[i]-y[i-1],2));
					NX[i] = N[0]; NY[i] = N[1];
				}

				double[] p = extract_profile(x[i], y[i], N);
				for (int j=0;j<w_profile;j++)
				{
					profile[i][j] = p[j];
				}
			}
        } 
		else
		{
			profile = new double[0][0];
		}
	return profile;
	}

	/******************************************************************************************
	 * Resample ROI data for refining the path
	 */
	private double[][] Refine_Profile(double[][] profile0)
	{
		// characterstics of the new profile
		float Ltot = L[L.length-1];
		int ntot = (int)(Ltot/l_profile + 0.5);
		int index = 0;
		int lcurrent = 0;
		if (ntot > 2.*L.length)
		{
			// New data structure
			double[][] profile = new double[ntot+1][w_profile];
			float[] LL = new float[ntot+1];
			float[] X = new float[ntot+1];
			float[] Y = new float[ntot+1];
			float[] VNX = new float[ntot+1];
			float[] VNY = new float[ntot+1];
					
			// Initial conditions
			LL[0] = 0;
			X[0] = x[0];
			Y[0] = y[0];
	
			// Get new profile
			float[] T0 = get_tangent(x[index],y[index],x[index+1],y[index+1]);
			float[] N0 ={0,0}; N0[0] = -T0[1];  N0[1] = T0[0];
			
			double[] p0 = extract_profile(X[0], Y[0], N0);
			for (int j=0;j<w_profile;j++)
			{
				profile[0][j] = p0[j];
				}		
				
			for (int i=1;i<=ntot;i++)
			{
				LL[i] = LL[i-1] + l_profile;
				if (index < x.length-1)
				{
					if(LL[i] > L[index+1]) // Start of a new coarse segment
					{
						index +=1;
						if (index == x.length-1)			// Last point of the path
						{
							float[] T = get_tangent(x[index-1],y[index-1],x[index],y[index]);
							float[] N ={0,0}; N[0] = -T[1];  N[1] = T[0];
							float dl = LL[i] - L[index];
	
							// Get new coordinates
							X[i] = 	x[index];
							Y[i] = 	y[index];
	
							// Get new profile
							double[] p = extract_profile(X[i], Y[i], N);
							for (int j=0;j<w_profile;j++)
							{
								profile[i][j] = p[j];
							}
							//break;						
						}
						else			// Entering a new segment of the original path
						{	
							float[] T = get_tangent(x[index],y[index],x[index+1],y[index+1]);
							float[] N ={0,0}; N[0] = -T[1];  N[1] = T[0];
							float dl = LL[i] - L[index];
							// Get new coordinates
							X[i] = 	x[index] + dl*T[0];
							Y[i] = 	y[index] + dl*T[1];
							
							// Get new profile
							double[] p = extract_profile(X[i], Y[i], N);
							for (int j=0;j<w_profile;j++)
							{
								profile[i][j] = p[j];
							}
						}
					}
					else			// keep going on an existing coarse segment of the original path
					{
						float[] T = get_tangent(x[index],y[index],x[index+1],y[index+1]);
						float[] N ={0,0}; N[0] = -T[1];  N[1] = T[0];
						float dl = LL[i] - L[index];
						
						// Get new coordinates
						X[i] = 	x[index] + dl*T[0];
						Y[i] = 	y[index] + dl*T[1];
						// Get new profile
						double[] p = extract_profile(X[i], Y[i], N);
						
						for (int j=0;j<w_profile;j++)
						{
							profile[i][j] = p[j];
						}
					}	 
				}
				
			}
		// Add the last point
		float[] T = get_tangent(x[x.length-2],y[x.length-2],x[x.length-1],y[x.length-1]);
		float[] N ={0,0}; N[0] = -T[1];  N[1] = T[0];
		X[X.length-1] = x[x.length-1];
		Y[X.length-1] = y[x.length-1];
		double[] p = extract_profile(X[X.length-1], Y[X.length-1], N);
		for (int j=0;j<w_profile;j++)
		{
			profile[X.length-2][j] = p[j];
		}
	
		// update
		L = LL;
		x = X;
		y = Y;
		return profile;
	}
	else
	{
		return	 profile0;
	}
	}
	/******************************************************************************************
	 * Export Data
	 */
	private void Refresh_Roi()
	{
		imp.deleteRoi();
		int[] XXi = new int[x.length];
		int[] YYi = new int[y.length];
		for (int i=0;i<XXi.length;i++)
		{
			XXi[i] = (int)(x[i]+0.5);
			YYi[i] = (int)(y[i]+0.5);
		}
		PolygonRoi Proi = new PolygonRoi(XXi,YYi,XXi.length,Roi.POLYLINE);
		imp.setRoi((Roi)Proi);
	}
	/******************************************************************************************
	 * Export Data
	 */
	private void Export_Data (double[][] profile)
	{
		if (profile.length >2)
		{
			String header = "Dist\t";
			for (int i=0;i<w_profile;i++) {header = header + "I" + i + "\t"; }
			IJ.log(header);
			for (int i=0;i<profile.length;i++)
			{
				String line = "" + L[i]+ "\t";
				for (int j=0;j<w_profile;j++)
				{
					line = line + profile[i][j]+ "\t";
				}
				IJ.log(line);
			}
		}		
	}

	/******************************************************************************************
	 * Extract profile along the axis perpendicular to the root
	 */
	private double[] extract_profile(float xx, float yy, float[] N)
		{
		double[] profile = new double[w_profile];
		for (int j=0;j<w_profile;j++)
		{
			int xj = (int)xx+(int)((j-w_profile/2.)*N[0]+0.5);
			int yj = (int)yy+(int)((j-w_profile/2.)*N[1] +0.5);
			profile[j] = 0;
			if (xj>=0 & yj>=0 & xj<width & yj<height)
				{profile[j] = I[xj][yj];}
		}
		return profile;
		}
	/******************************************************************************************
	 * Get unit tangent vector
	 */	
	private float[] get_tangent(float x1, float y1, float x2, float y2)
		{
		float[] TT = new float[2];
		TT[0] = x2 - x1;
		TT[1] = y2 - y1;
		float l = (float)Math.sqrt(TT[0]*TT[0] + TT[1]*TT[1]);
		TT[0] = TT[0]/l;
		TT[1] = TT[1]/l;
	
		return TT;
	}

}


