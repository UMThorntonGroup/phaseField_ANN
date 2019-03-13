/*

ic_generator.i

//Works for periodic boundary conditions ONLY
generates an array of random positions and radii of nuclei in 2D or 3D
in the following format
Positions:
 {{0.1,0.3,0},{0.8,0.7,0},{0.5,0.2,0},{0.4,0.4,0},{0.3,0.9,0},...,{1,1,0},{0.7,0.95,0}};
 Radii:
 {12, 14, 19, 16, 11, 12, 17, ... , 11, 14};
 
*/
#include "random.i"

//MAIN
//Defining constants

//dimensionality
dim=3;

//System size
lx=512.0;
ly=512.0;
lz=512.0;

lvec=[lx,ly,lx];

//Number of seeds
nseeds = 120;
//Average radius of seeds
nrad = 5.0;
//Std. deviation of seed radii around average
sigma = 0.5;
//Minimum distance between nuclei
mindist = 3.0;

lvec=lvec(1:dim);

random_seed, 0.12345678;

//Initializing position and dimension vectors;
nucposvec=array(0.0,[2,dim,nseeds]);
nucsizevec=array(0.0,nseeds);

//First nucleus
nucpos = lvec*random(dim);
nucsize = nrad + sigma*random_n();
write, "Nucleus: ", 1;
write, "Nucleus at ";
nucpos;
write, " with radius ", nucsize;
nucposvec(,1)=nucpos;
nucsizevec(1)=nucsize;

for (j=2;j<=nseeds;j++){
    accept=0;
    for (jc=1; accept==0 ;jc++){
        nucpos = lvec*random(dim);
        nucsize = nrad + sigma*random_n();
        write, "Nucleus: ", j;
        write, "Trying new nucleus at ";
        nucpos;
        write, " with radius ", nucsize;
        //Check no overlap against previous nucleus
        reject=0;
        for (k=1;k<=j-1;k++){
            //Minimum distance for periodic boundary conditions
            shortestdv=array(0.0,dim);
            for (jd=1;jd<=dim;jd++){
                distcomp=abs(nucpos(jd)-nucposvec(jd,k));
                shortestdv(jd)=min(distcomp,lvec(dim)-distcomp);
            }
            dist_to_new_nuc = sqrt(sum(shortestdv^2.0))-nucsize-nucsizevec(k);
            //write, "Distance to nucleus ", k, " = ", dist_to_new_nuc;
            if (dist_to_new_nuc < mindist) reject=1;
        }
        if (reject==1) write, "nucleus rejected";
        if (reject==0) accept=1;
    }
nucposvec(,j)=nucpos;
nucsizevec(j)=nucsize;
}

if (dim==2){
    winkill, 0;
    window, 0, dpi=100;
    limits, 0, lx, 0, ly;
    fma;
    for(j=1;j<=nseeds;j++){
        plmk, nucposvec(2,j),nucposvec(1,j), marker=4, width=10;
        //rdline;
    }
}

//Writing out positions readable in icc
fl=create("nuclei_data.txt");
write, fl, "Nuclei coordinates";
write, fl, format="%s", "{";
for(j=1;j<=nseeds;j++){
    if (j<nseeds){
        if (j%4) {
            write, fl, format="%s", "{";
            for (k=1;k<=dim;k++){
                if (k<dim) write, fl, format="%.2f, ", nucposvec(k,j);
                if (k==dim) write, fl, format="%.2f", nucposvec(k,j);
            }
            write, fl, format="%s", "},";
        } else {
            write, fl, format="%s", "{";
            for (k=1;k<=dim;k++){
                if (k<dim) write, fl, format="%.2f, ", nucposvec(k,j);
                if (k==dim) write, fl, format="%.2f", nucposvec(k,j);
            }
            write, fl, format="%s\n", "},";
        }
    }
    if (j==nseeds){
        write, fl, format="%s", "{";
        for (k=1;k<=dim;k++){
            if (k<dim) write, fl, format="%.2f, ", nucposvec(k,j);
            if (k==dim) write, fl, format="%.2f", nucposvec(k,j);
        }
        write, fl, format="%s", "}";
    }
}
write, fl, format="%s \n\n", "};";

write, fl, "Nuclei radii";
write, fl, format="%s", "{";
for(j=1;j<=nseeds;j++){
    if (j<nseeds){
        if (j%4) {
            write, fl, format="%.3f, ", nucsizevec(j);
        } else {
            write, fl, format="%.3f, \n", nucsizevec(j);
        }
    }
    if (j==nseeds) write, fl, format="%.3f", nucsizevec(j);
}
write, fl, format="%s", "};";
close, fl;

