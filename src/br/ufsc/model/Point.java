/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package br.ufsc.model;

import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;

/**
 *
 * @author vanes
 */
public class Point {

    protected MultipleAspectTrajectory trajectory;
    protected int rId;
    protected double x;
    protected double y;
    protected List<AttributeValue> listAttrValues;
    protected TemporalAspect time;
    protected String cellReference;
//    public static SimpleDateFormat formatDate = new SimpleDateFormat("HH:mm");
    
    

    public Point(MultipleAspectTrajectory t, int rid, double x, double y, Date startTime) {
        this.trajectory = t;
        this.rId = rid;
        this.x = x;
        this.y = y;
        this.time = new TemporalAspect(startTime);
        listAttrValues = new ArrayList<>();
    }
    
        public Point(MultipleAspectTrajectory t, int rid, double x, double y, Date startTime, Date endTime) {
        this.trajectory = t;
        this.rId = rid;
        this.x = x;
        this.y = y;
        this.time = new TemporalAspect(startTime,endTime);
        listAttrValues = new ArrayList<>();
    }

    public Point(int rid, double x, double y, Date startTime) {
        this.rId = rid;
        this.x = x;
        this.y = y;
        listAttrValues = new ArrayList<>();
        this.time = new TemporalAspect(startTime);
    }
    
        public Point(int rid, double x, double y, Date startTime, Date endTime) {
        this.rId = rid;
        this.x = x;
        this.y = y;
        listAttrValues = new ArrayList<>();
        this.time = new TemporalAspect(startTime, endTime);
    }

    public Point(double x, double y, Date startTime) {
        this.x = x;
        this.y = y;
        this.time = new TemporalAspect(startTime);
        listAttrValues = new ArrayList<>();
    }
    
    public Point(double x, double y, Date startTime, Date endTime) {
        this.x = x;
        this.y = y;
        this.time = new TemporalAspect(startTime, endTime);
        listAttrValues = new ArrayList<>();
    }
    

    public Point(double x, double y, int startTime) {
        this.x = x;
        this.y = y;
        this.time = new TemporalAspect(Util.convertMinutesToDate(startTime));
        listAttrValues = new ArrayList<>();
    }
    
        public Point(double x, double y, int startTime, int endTime) {
        this.x = x;
        this.y = y;
        this.time = new TemporalAspect(Util.convertMinutesToDate(startTime), Util.convertMinutesToDate(endTime));
        listAttrValues = new ArrayList<>();
    }

    public Point(int rid, double x, double y, int startTime) {
        this.rId = rid;
        this.x = x;
        this.y = y;
        this.time = new TemporalAspect(Util.convertMinutesToDate(startTime));
        listAttrValues = new ArrayList<>();
    }
    
        public Point(int rid, double x, double y, int startTime, int endTime) {
        this.rId = rid;
        this.x = x;
        this.y = y;
        this.time = new TemporalAspect(Util.convertMinutesToDate(startTime), Util.convertMinutesToDate(endTime));
        listAttrValues = new ArrayList<>();
    }

    public Point(int rid, double x, double y, int startTime, ArrayList semantics) {
        this.rId = rid;
        this.x = x;
        this.y = y;
        this.time = new TemporalAspect(Util.convertMinutesToDate(startTime));
        listAttrValues = semantics;
    }
    
    public Point(int rid, double x, double y, int startTime, int endTime, ArrayList semantics) {
        this.rId = rid;
        this.x = x;
        this.y = y;
        this.time = new TemporalAspect(Util.convertMinutesToDate(startTime), Util.convertMinutesToDate(endTime));
        listAttrValues = semantics;
    }
    
    public Point(int rid, double x, double y, Date startTime, ArrayList semantics) {
        this.rId = rid;
        this.x = x;
        this.y = y;
        this.time = new TemporalAspect(startTime);
        listAttrValues = semantics;
    }
    public Point(int rid, double x, double y, Date startTime, Date endTime, ArrayList semantics) {
        this.rId = rid;
        this.x = x;
        this.y = y;
        this.time = new TemporalAspect(startTime, endTime);
        listAttrValues = semantics;
    }
    
    
    
    public Point(double x, double y) {
        this.x = x;
        this.y = y;
        listAttrValues = new ArrayList<>();
    }
    
    public Point(){
        listAttrValues = new ArrayList<>();
    }
    
    public double euclideanDistance(Point other) {
        return Math.sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
    }

    @Override
    public String toString() {
        Locale.setDefault(Locale.US);
        DecimalFormat formatNumber = new DecimalFormat("##.##");
        String txt = ((rId != 0) ? ("RId: " + rId + "\n") : ("rt"));

        txt += "\n(x,y)= (" + formatNumber.format(x) + "," + formatNumber.format(y) + "), ";
        if (time != null) {
//            txt += "\nTIME: " + formatDate.format(time) + ", ";
            txt += "\nTemporal Aspect: " + time + ", ";
        }
        if (!listAttrValues.isEmpty() ) {
            txt += showAttrValues();
        }
//                );
        return txt;

    }

    public MultipleAspectTrajectory getTrajectory() {
        return trajectory;
    }

    public void setTrajectory(MultipleAspectTrajectory trajectory) {
        this.trajectory = trajectory;
    }

    public int getrId() {
        return rId;
    }

    public void setrId(int rId) {
        this.rId = rId;
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public double getY() {
        return y;
    }

    public void setY(double y) {
        this.y = y;
    }

    public TemporalAspect getTime() {
        return time;
    }

    public void setTime(Date startTime) {
        this.time = new TemporalAspect(startTime);
    }
    
    public void setTime(Date startTime, Date endTime) {
        this.time = new TemporalAspect(startTime, endTime);
    }
    
    
        
    
    public void addAttrValue(Object value, SemanticAspect attr) {
        listAttrValues.add(new AttributeValue(value, attr));
    }
    
        
    public void addAttrValue(Object value, SemanticAspect attr, double numValueSD) {
        listAttrValues.add(new AttributeValue(value, attr, numValueSD));
    }

    public String showAttrValues() {
        String txt = "(";
        for (AttributeValue atv : listAttrValues) {
                 
            txt += atv.getAttibute().getName() + ": " + atv.getValue() + ", ";
        }
        txt += ")";

        return txt;
    }
    
    

    public List<AttributeValue> getListAttrValues() {
        return listAttrValues;
    }

    public void setListAttrValues(List<AttributeValue> listAttrValues) {
        this.listAttrValues = listAttrValues;
    }

    public int getTimeInMinutes(Date t){
        Calendar c = Calendar.getInstance();
        c.setTime(t);
        return (c.get(Calendar.HOUR_OF_DAY) * 60) + c.get(Calendar.MINUTE);
    }
    
    public int getTimeInMinutes(){
        Calendar c = Calendar.getInstance();
        c.setTime(time.getStartTime());
        return (c.get(Calendar.HOUR_OF_DAY) * 60) + c.get(Calendar.MINUTE);
    }
    
    /**
     * Return the value of the attribute a
     * @param SemanticAspect attribute
     * @return value of the attribute
     */
    public AttributeValue getAttributeValue(SemanticAspect attribute){
        for(AttributeValue atv: listAttrValues){
            if(atv.getAttibute().equals(attribute))
                return atv;
        }
        return null;
    }
    
    
    // Equals and HashCode

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 89 * hash + Objects.hashCode(this.trajectory);
        hash = 89 * hash + this.rId;
        hash = 89 * hash + (int) (Double.doubleToLongBits(this.x) ^ (Double.doubleToLongBits(this.x) >>> 32));
        hash = 89 * hash + (int) (Double.doubleToLongBits(this.y) ^ (Double.doubleToLongBits(this.y) >>> 32));
        hash = 89 * hash + Objects.hashCode(this.listAttrValues);
        hash = 89 * hash + Objects.hashCode(this.time);
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
        final Point other = (Point) obj;
        if (this.rId != other.rId) {
            return false;
        }
        if (Double.doubleToLongBits(this.x) != Double.doubleToLongBits(other.x)) {
            return false;
        }
        if (Double.doubleToLongBits(this.y) != Double.doubleToLongBits(other.y)) {
            return false;
        }
        if (!Objects.equals(this.trajectory, other.trajectory)) {
            return false;
        }
        if (!Objects.equals(this.listAttrValues, other.listAttrValues)) {
            return false;
        }
        if (!Objects.equals(this.time, other.time)) {
            return false;
        }
        return true;
    }

    public String getCellReference() {
        return cellReference;
    }

    public void setCellReference(String cellReference) {
        this.cellReference = cellReference;
    }
    
    
    
}
