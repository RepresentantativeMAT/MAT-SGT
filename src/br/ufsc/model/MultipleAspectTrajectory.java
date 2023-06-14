/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package br.ufsc.model;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author vanes
 */
public class MultipleAspectTrajectory implements Cloneable {

    

    private String description;
    private int id;
    private List<Point> pointList;
    private int coverPoints;
    private boolean dailyInfo;

    public MultipleAspectTrajectory(String description, int id) {
        this.description = description;
        this.id = id;
        pointList = new ArrayList<Point>();
    }

    public MultipleAspectTrajectory(int id) {
        this.id = id;
        pointList = new ArrayList<Point>();
    }

    public MultipleAspectTrajectory(String description) {
        this.description = description;
        pointList = new ArrayList<Point>();
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public void addPoint(Point p) {
        pointList.add(p);
        p.setTrajectory(this);
    }

    public void removePoint(Point p) {
        pointList.remove(p);
        p.setTrajectory(null);
    }

    public Point getLastPoint() {
        return pointList.get(pointList.size() - 1);
    }

    public int getCoverPoints() {
        return coverPoints;
    }

    public void setCoverPoints(int coverPoints) {
        this.coverPoints = coverPoints;
    }

    /**
     * Method to increment some data points (size of points) to compute coverPoints,
     * -- i.e., the input data points covered by our representative trajectory
     * @param sizeDataPoints 
     */
    public void incrementValue(int sizeDataPoints) {
//        System.out.println("Cover points: " + coverPoints + " plus: " + sizeDataPoints);
        this.coverPoints += sizeDataPoints;
    }
    
    public void decrementValue(int value) {
//        System.out.println("Cover points: " + coverPoints + " minus: " + value);
        this.coverPoints -= value;
    }

    @Override
    public String toString() {
        String aux = "ID: " + id;
        aux += "\nDescription: " + description;
        if (!pointList.isEmpty()) {
            aux += "\nPoint List: \n";

            for (Point p : pointList) {
                aux += p + "\n";
            }
        }
        return aux;
    }

    public List<Point> getPointList() {
        return pointList;
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    public boolean isDailyInfo() {
        return dailyInfo;
    }

    public void setDailyInfo(boolean dailyInfo) {
        this.dailyInfo = dailyInfo;
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 97 * hash + this.id;
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final MultipleAspectTrajectory other = (MultipleAspectTrajectory) obj;
        return this.id == other.id;
    }
}
