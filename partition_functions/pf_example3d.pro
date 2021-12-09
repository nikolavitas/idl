PRO pf_example3d

z = 25

nx = 256
ny = 128
nz = 256
a = FLTARR(nx, ny, nz)
OPENR, 1, '~/linux/chemeq_v2.0/t3d.003003.dat'
READU, 1, a
CLOSE, 1
t = a

;pfi = PF(t, z, data = 'i')
;pfw = PF(t, z, data = 'w')
;pfg = PF(t, z, data = 'g')
;pfs = PF(t, z, data = 's')
pfc = PF(t, z, data = 'c')

PRINT, 'Gray:'
PRINT, pfg.u1
PRINT, pfg.u2
PRINT, pfg.u3

PRINT, 'Wittmann:'
PRINT, pfw.u1
PRINT, pfw.u2
PRINT, pfw.u3

PRINT, 'Irwin:'
PRINT, pfi.u1
PRINT, pfi.u2
PRINT, pfi.u3

PRINT, 'Cowley:'
PRINT, pfc.u1
PRINT, pfc.u2
PRINT, pfc.u3
PRINT, pfc.u4

PRINT, 'Sauval and Tatum:'
PRINT, pfs.u1
PRINT, pfs.u2

!p.background = GETCOLOR('white')
!p.color = GETCOLOR('black')
WINDOW, xsi = 500, ysi = 400
PLOT,  t,pfg.u1, lines = 0, thick = 2, title = "Partition function for neutral Mn", $
       pos = [0.10, 0.14, 0.94, 0.92], xtit = 'T [K]', ytit = 'U', charsi = 1.4
OPLOT, t,pfw.u1, lines = 0, color = getcolor('Orange Red'), thick = 2
OPLOT, t,pfi.u1, lines = 0, color = getcolor('Steel Blue'), thick = 2
OPLOT, t,pfc.u1, lines = 0, color = getcolor('Forest Green'), thick = 2
OPLOT, t,pfs.u1, lines = 0, color = getcolor('Chocolate'), thick = 2
OPLOT, [2400, 3400], [35, 35], thick = 2
OPLOT, [2400, 3400], [33, 33], color = getcolor('Orange Red'), thick = 2
OPLOT, [2400, 3400], [31, 31], color = getcolor('Steel Blue'), thick = 2
OPLOT, [2400, 3400], [29, 29], color = getcolor('Forest Green'), thick = 2
OPLOT, [2400, 3400], [27, 27], color = getcolor('Chocolate'), thick = 2
XYOUTS, 3450, 34.4, 'Gray', chars = 1.4
XYOUTS, 3450, 32.4, 'Wittmann', chars = 1.4
XYOUTS, 3450, 30.4, 'Irwin', chars = 1.4
XYOUTS, 3450, 28.4, 'Cowley', chars = 1.4
XYOUTS, 3450, 26.4, 'Sauval and Tatum', chars = 1.4
WRITE_PNG, 'example_atomic_pf_mn.png', tvrd(/true)



OPLOT, t, pfg.u2, lines = 1
OPLOT, t, pfg.u3, lines = 2

OPLOT, t,pfw.u2, lines = 1, color = getcolor('red')
OPLOT, t,pfw.u3, lines = 2, color = getcolor('red')

OPLOT, t,pfi.u2, lines = 1, color = getcolor('blue')
OPLOT, t,pfi.u3, lines = 2, color = getcolor('blue')

OPLOT, t,pfc.u2, lines = 1, color = getcolor('green')
OPLOT, t,pfc.u3, lines = 2, color = getcolor('green')

WRITE_PNG, 'example_abundances_2.png', tvrd(/true)


stop
END
