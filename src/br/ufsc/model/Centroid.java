/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package br.ufsc.model;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;

/**
 *
 * @author vanes
 */
public class Centroid extends Point implements Comparable{

    private List<Point> pointListSource;
    private List<STI> listSTI;
    private STI sti;
//    Locale.setDefault(Locale.US);
//    DecimalFormat formatNumber = new DecimalFormat("0.###");

    public Centroid() {
        super();
        pointListSource = new ArrayList<>();
        listSTI = new ArrayList<>();
    }

    public Centroid(double x, double y) {
        super(x, y);
        pointListSource = new ArrayList<>();
        listSTI = new ArrayList<>();
    }

    public Centroid(Point p) {
        super(p.getX(), p.getY());
        pointListSource = new ArrayList<>();
        listSTI = new ArrayList<>();
    }

    //methods
    
    public void setSpatialDimension(double x, double y) {
        this.x = x;
        this.y = y;
    }

    
    public STI getSti() {
        return sti;
    }

    public void setSti(STI sti) {
        this.sti = sti;
        if(sti != null)
            this.sti.setPoint(this);
    }

    
    
    public void addPoint(Point p) {
        pointListSource.add(p);
    }

    public void removePoint(Point p) {
        pointListSource.remove(p);
    }

    public void addSTI(STI sti) {
        sti.setPoint(this);
        listSTI.add(sti);
        
    }

    public void addSTI(Date startTime, float proportion) {
//        sti.setPoint(this);
        listSTI.add(new STI(startTime, proportion, this));
    }

    public void addSTI(Date startTime, Date endTime, float proportion) {
//        sti.setPoint(this);
        listSTI.add(new STI(startTime, endTime, proportion, this));
    }

    public void removeSTI(STI time) {
        time.setPoint(null);
        listSTI.remove(time);
    }
    
    
    

//    @Override
    public String toString() {
        Locale.setDefault(Locale.US);
        DecimalFormat formatNumber = new DecimalFormat("##.##");
        String aux = ("rt");

        aux += "\n(x,y)= (" + formatNumber.format(x) + "," + formatNumber.format(y) + "), ";
//        if(!listSTI.isEmpty()){
//            aux += "\nTime: "+listSTI;
//        }
        if (!listAttrValues.isEmpty() ) {
            aux += showAttrValues();
        }
        if(getSti()!=null){
            aux += "\nTime: "+sti;
        }
                aux += "\nCell: "+cellReference
                        +"\nLocal mapped: ";
        for (Point p : pointListSource) {
            aux += p.getTrajectory().getId() + " - " + p.getrId() + ", ";
        }
        aux += "";
        return aux;
    }

    public List<Point> getPointListSource() {
        return pointListSource;
    }

    @Override
    public String showAttrValues() {

        Locale.setDefault(Locale.US);
        DecimalFormat formatNumber = new DecimalFormat("0.###");

        String txt = "";
        for (AttributeValue atv : super.listAttrValues) {
            txt += "\n";

            if (!(atv.getValue() instanceof Map)) {
                txt += atv.getAttibute().getName() + "= " + atv.getValue() + "";
            } else {

                txt += "Ranking of " + atv.getAttibute().getName() + "= [";
                HashMap<String, Double> allValues = (HashMap) atv.getValue();

                for (Map.Entry<String, Double> eachValue : allValues.entrySet()) {
                    txt += eachValue.getKey().replace(",", ";") + " -> " + formatNumber.format(eachValue.getValue()) + ", ";
                }

                txt += " ], ";
            }
        }
//        System.out.println("Lista STI: "+listSTI);
        //txt += "";

        if (!listSTI.isEmpty()) {
            txt += "Ranking of Temporal Interval = [";
            for (STI sti : listSTI) {
                txt += sti + ", ";
            }
            txt += " ], ";
        }

        return txt;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 17 * hash + Objects.hashCode(this.pointListSource);
        hash = 17 * hash + Objects.hashCode(this.listSTI);
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
        final Centroid other = (Centroid) obj;
        
        if (!Objects.equals(this.listSTI, other.listSTI)) {
            return false;
        }
        return true;
    }

    
    @Override
    public int compareTo(Object other) {
        
        if(this.getSti().getInterval().getStartTime().after(((Centroid)other).getSti().getInterval().getStartTime()))
            return 1;
        else if (this.getSti().getInterval().getStartTime().before(((Centroid)other).getSti().getInterval().getStartTime()))
            return -1;
        else
            return 0;
        
    }


    
    

}
