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
public class Cytodex_Actin {

    // index
    private int index;
    // volume
    private double volume;
    // distances center object from global center
    private double centerDistToCenter;
    // distances border object to global center
    private double borderDistToCenter;
    // nucleus number
    private int nucNumber;

    
    
    
    public Cytodex_Actin(int index, double volume, double centerDistToCenter, double borderDistToCenter, int nucNumber) {
            this.index = index;
            this.volume = volume;
            this.centerDistToCenter = centerDistToCenter;
            this.borderDistToCenter = borderDistToCenter;
            this.nucNumber = nucNumber;
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
    
    public void setBorderDistToCenter(double borderDistToCenter) {
        this.borderDistToCenter = borderDistToCenter;
    }
    
    public double getBorderDistToCenter() {
        return borderDistToCenter;
    }
    
     public void setNucNumber(int nucNumber) {
        this.nucNumber = nucNumber;
    }
    
    public int getNucNumber() {
        return nucNumber;
    }
    
    
}
