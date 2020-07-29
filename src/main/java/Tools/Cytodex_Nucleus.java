/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Tools;

/**
 *
 * @author phm
 */
public class Cytodex_Nucleus {
    
    
     // index
    private int index;
    // volume
    private double volume;
    // distances center object from global center
    private double centerDistToCenter;
    // Sphericity
    private double sphericity;
    // Compactness
    private double compactness;
    // angle between centroid and main axes
    private double centroidAngle;
    // angle between main axes and all nucleus main axe
    private double allAngle;
    // nucleus number
    private int actinObjectIndex;
    // distance to closest nucleus
    private double nucClosestDist;
    // actin diameter
    private double actinDiameter;
    
    
     public Cytodex_Nucleus(int index, double volume, double centerDistToCenter, double sphericity, double compactness, double centroidAngle, double allAngle,
             int actinObjectIndex, double nucClosestDist, double actinDiameter) {
            this.index = index;
            this.volume = volume;
            this.centerDistToCenter = centerDistToCenter;
            this.sphericity = sphericity;
            this.compactness = compactness;
            this.centroidAngle = centroidAngle;
            this.allAngle = allAngle;
            this.actinObjectIndex = actinObjectIndex;
            this.nucClosestDist = nucClosestDist;
            this.actinDiameter = actinDiameter;
    }
    
    public void setIndex(int index) {
        this.index = index;
    }
    
    public int getIndex() {
        return index;
    }
    
    public void setVolume(double volume) {
        this.volume = volume;
    }
    
    public double getVolume() {
        return volume;
    }

    public void setCenterDistToCenter(double centerDistToCenter) {
        this.centerDistToCenter = centerDistToCenter;
    }
    
    public double getCenterDistToCenter() {
        return centerDistToCenter;
    }
    
    public void setNucClosestDist(double centerDistToCenter) {
        this.nucClosestDist = nucClosestDist;
    }
    
    public double getNucClosestDist() {
        return nucClosestDist;
    }
    
    public void setSphericity(double sphericity) {
        this.sphericity = sphericity;
    }
    
    public double getSphericity() {
        return sphericity;
    }
    
    public void setCompactness(double compactness) {
        this.compactness = compactness;
    }
    
    public double getCompactness() {
        return compactness;
    }
    
    public void setCentroidAngle(double centroidAngle) {
        this.centroidAngle = centroidAngle;
    }
    
    public double getCentroidAngle() {
        return centroidAngle;
    }
    
    public void setAllAngle(double allAngle) {
        this.allAngle = allAngle;
    }
    
    public double getAllAngle() {
        return allAngle;
    }
    public void setActinObjectIndex(int actinObjectIndex) {
        this.actinObjectIndex = actinObjectIndex;
    }
     
    public int getActinObjectIndex() {
        return actinObjectIndex;
    }
    
    public void setActinDiameter(double actinDiameter) {
        this.actinDiameter = actinDiameter;
    }
    
    public double getActinDiameter() {
        return actinDiameter;
    }
}
