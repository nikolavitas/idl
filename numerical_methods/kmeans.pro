;+
; NAME:
;
;       KMEANS
;
; PURPOSE:
;
;       This function performs k-means clustering of a given dataset for a
;       given number of clusters. Three initializations are available:
;       completely random (standard k-means), the furthest point, and the 
;       k-means++ (Arthur & Vassilvitskii (2007), k-means++ : The Advantages 
;       of Careful Seeding).
; 
;       The function aslo computes silhouette of the final clustering, see 
;       Rousseeuw (1987), "Silhouettes: a Graphical Aid to the Interpretation 
;       and Validation of Cluster Analysis". The silhouette is computed for 
;       every measurement. The value averaged over the entire data set is a 
;       measure how appropriate is the choice of the clusters. See more at 
;       https://en.wikipedia.org/wiki/Silhouette_(clustering)
; 
;
;       IDL has its own native routine for k-means clustering: CLUSTER. Note
;       that an earlier version of that routine was called KMEANS, but became
;       obsolete. My routine is not vectorized, but it is more readable and 
;       thus easier to change and upgrade. It also has features missing in 
;       the native CLUSTER function.
;      
; AUTHOR:
;
;       Nikola Vitas
;       Instituto de Astrofisica de Canarias (IAC)
;       C/ Via Lactea, s/n
;       E38205 - La Laguna (Tenerife), Espana
;       Email: n.vitas@iac.es, nikola.vitas@gmail.com
;       Homepage: nikolavitas.blogspot.com
;
; CATEGORY:
;
;       Data mining, statistics
;
; CALLING SEQUENCE:
;
;       y = KMEANS(x, k = k, initial = initial)
;
; INPUTS:
;
;       z       = Float 2D array. All data in the matrix, m x n, where m is
;                 the dimension and n is the number of measurements.
; 
; OUTPUTS:
;
;       y       = Structure with tags: clusters (number of assigned cluster for 
;                 each measurement), frequencies (percentual fraction of the 
;                 total number of measurements that is assigned to each of the 
;                 clusters), and silhouette (measure of how appropriate is the 
;                 final clustering for the given data)
;  
; INPUT KEYWORDS:
;
;       k       = Number of clusters.
; 
;       initial = String. Initialization scheme: ['random'|'furthest'|'kmeanspp']
;
; OUTPUT KEYWORDS:
;
; EXAMPLE:
;
;       IDL> n = 100   ; Number of points
;       IDL> m = 2     ; Number of dimensions, ie. length of vectors
;       IDL> k = 4     ; Number of clusters
;       IDL>
;       IDL> x = REFORM(RANDOMU(10, m*n)/4., [m, n])
;       IDL> x[0, 0:24] += 0.2
;       IDL> x[1, 0:24] += 0.2
;       IDL> x[0, 25:49] += 0.8
;       IDL> x[1, 25:49] += 0.8
;       IDL> x[0, 50:74] += 0.2
;       IDL> x[1, 50:74] += 0.8
;       IDL> x[0, 75:99] += 0.8
;       IDL> x[1, 75:99] += 0.2
;       IDL> y = KMEANS(x, k = k, initial = 'kmeanspp') 
;
;       ; Show the initial data points
;       IDL> PLOT, x[0, *], x[1, *], psym = 4
;
;       ; Show which are assigned to the cluster 2
;       IDL> id = WHERE(y.clusters EQ 2)
;       IDL> OPLOT, x[0, id], x[1, id], psym = -1
;      
; DEPENDENCIES:
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, March 2017
;-
;================================================================================
; KMEANS, IDL routine by Nikola Vitas is licensed under a Creative Commons 
; Attribution-NonCommercial-ShareAlike 3.0 Unported License.
;
; This software is provided by NV ''as is'' and any express or implied warranties, 
; including, but not limited to, the implied warranties of merchantability and 
; fitness for a particular purpose are disclaimed. In no event shall NV be liable 
; for any direct, indirect, incidental, special, exemplary, or consequential 
; damages (including, but not limited to, procurement of substitute goods or 
; services; loss of use, data, or profits; loss of use, data, or profits; or 
; business interruption) however caused and on any theory of liability, whether 
; in contract, strict liability, or tort (including negligence or otherwise) 
; arising in any way out of the use of this software, even if advised of the 
; possibility of such damage.
;================================================================================

FUNCTION kmeans, x, k = k, initial = initial

seed = 7
m = (SIZE(x))[1] ; number of dimensions 
n = (SIZE(x))[2] ; number of measurements

clusters0 = INTARR(n)  

; Cluster centers
r = FLTARR(m, k)  

IF initial EQ 'random' THEN BEGIN
  ; Set the random seeds:
  idseed = LONG((RANDOMU(seed, k))*LONG(n))
  
  FOR l = 0, k-1 DO r[*, l] = x[*, idseed[l]]
ENDIF

IF initial EQ 'furthest' THEN BEGIN

  ; Set the first random seed:
  idseed = LONG((RANDOMU(seed, 1))*LONG(n))

  ; Set the first random seed:
  r[*, 0] = x[*, idseed]
  
  ; Compute distances of all data points from each seed and find the furthest 
  ; one.
  FOR l = 1, k-1 DO BEGIN
    d = FLTARR(n)
    FOR l1 = 1, l DO $
      FOR j = 0L, n-1 DO d[j] += SQRT(TOTAL((x[*, j] - r[*, l1-1])^2))
  
    id = WHERE(MAX(d) EQ d)
    r[*, l] = x[*, id]

  ENDFOR
ENDIF

IF initial EQ 'kmeanspp' THEN BEGIN

  ; Set the first random seed:
  idseed = LONG((RANDOMU(seed, 1))*LONG(n))

  ; First cluster center    
  r[*, 0] = x[*, idseed]
  
  d = FLTARR(k, n)
  
  FOR l = 0, k-2 DO BEGIN
  
    ; Distances from r[*, 0] to all data points
    FOR j = 0L, n-1 DO $ 
      d[l, j] = SQRT(TOTAL((x[*, j] - r[*, l])^2))
  
    dmin = MIN(d[0:l, *], dim = 1)
  
    ; Make probability distribution
    p = REFORM(dmin^2/TOTAL(dmin^2)) 
    
    ; Save the original order
    sorted = SORT(p)
    
    ; Sort probabilities
    p = p[sorted]  

    ; Compute cumulative probability
    f = TOTAL(p, /cum)
    
    ; Get random value with uniform distribution
    u = RANDOMU(seed, 1)
    
    ; Get random value with D^2 distribution
    xprim = INTERPOL(FINDGEN(n), f, u)
  
    ; Translate it into the original index.
    id = sorted[FIX(xprim)]
    
    ; Set the new cluster center
    r[*, l+1] = x[*, id]
        
  ENDFOR
    
ENDIF

nchanges = 100
niter = 0
WHILE nchanges GT 0 DO BEGIN    
  
  clusters = INTARR(n)
  d = FLTARR(k, n)
        
  ; Compute Euclidean distance between measurements and cluster centers
  FOR l = 0, k-1 DO $
    FOR j = 0L, n-1 DO BEGIN
      d[l, j] = SQRT(TOTAL((x[*, j] - r[*, l])^2))
    ENDFOR
  
  ; Update clasters
  FOR j = 0L, n-1 DO $
    clusters[j] = WHERE(d[*, j] EQ MIN(d[*, j]))  

  ; Update cluster centers  
  FOR l = 0, k-1 DO $
    FOR i = 0, m-1 DO $
      r[i, l] = TOTAL(x[i, *]*(clusters EQ l))/TOTAL(clusters EQ l)
  
  nchanges = total(clusters-clusters0 NE 0)
  clusters0 = clusters

  niter += 1
ENDWHILE

; Compute frequencies
frequencies = FLTARR(k)
FOR l = 0, k-1 DO frequencies[l] = TOTAL(clusters EQ l)/n*100

; Compute distance between any two measurements
d = FLTARR(n, n)
FOR j1 = 0L, n-1 DO $
  FOR j2 = 0L, n-1 DO $
    d[j1, j2] = SQRT(TOTAL((x[*, j1] - x[*, j2])^2))

; Compute dissimilarities and silhouette
a = FLTARR(n)
b = FLTARR(n)
s = FLTARR(n)
FOR j = 0L, n-1 DO BEGIN
  FOR l = 0, k-1 DO BEGIN
    IF clusters[j] EQ l THEN BEGIN
      a[j] = AVG(d[j, WHERE(clusters EQ l)])
      b[j] = AVG(d[j, WHERE(clusters NE l)])
      s[j] = (b[j] - a[j])/MAX([a[j], b[j]])
    ENDIF
  ENDFOR
ENDFOR
    
result = {clusters:clusters, frequencies:frequencies, silhouette:s}
    
RETURN, result 
END