!     ######################
!     #     tictoc.f90     #
!     ######################

!     Autor: Juan C. Toledo
!     Version: 1.2
!     Fecha: 3 Mar 2009

!     ############

!     Historial de version:
!     1.2: Se usa ahora SYSTEM_CLOCK en vez de DATE_AND_TIME para
!          calcular el tiempo transcurrido entre TIC y TOC
!     1.1: Se agrego STAMP para imprimir fecha y hora actuales
!     1.0: Version original

!     ############

!     Descripcion:

!     Contiene subrutinas que permiten medir el tiempo real
!     de ejecucion que transcurre entre dos puntos de un
!     programa.

!     ############

!     Contenido:

!     TIC(mark)

!       Marca el inicio del conteo.
!       El argumento mark es un INTEGER que se devuelve con el
!       tiempo transcurrido desde el inicio del programa
!       al momento de llamar la funcion. Sirve de punto
!       de referencia para llamadas subsecuentes a TOC.

!     TOC(mark)

!       Calcula e imprime el tiempo transcurrido.
!       El argumento mark es un INTEGER que contiene el tiempo
!       de un punto de ejecucion, previamente marcado por TIC.

!     STAMP()

!       Imprime la fecha y hora actuales

!     ############

!     Detalles:

!     La version 1.2 utiliza SYSTEM_CLOCK ya que tiene mas
!     resolucion y es mas simple de implementar que DATE_AND_TIME.
!     Sin embargo, SYSTEM_CLOCK solo puede contar hasta cierto
!     valor, que depende del compilador. Para el GCC 4.2.1, este
!     valor corresponde a 24.855 dias. Este programa usa el reloj
!     del sistema para medir el tiempo (wall time).

!     Para empezar a contar el tiempo se debe llamar la
!     subrutina TIC, dandole una variable INTEGER como
!     argumento. Se pueden usar multiples variables para
!     marcar distintos puntos en la ejecucion de un programa.

!     Para calcular el tiempo que ha transcurrido dese que se
!     marco cierto punto, se llama a la subrutina TOC dandole
!     como argumento la variable que se uso en TIC. La resolucion
!     puede variar de sistema a sistema pero en general es del
!     orden de milisegundos.

!     Para usar este archivo en un programa, es necesario incluirlo
!     añadiendo la siguiente linea al principio del codigo:

!     INCLUDE 'tictoc.f'

!     La subrutina STAMP imprime la fecha y hora actuales.

!     ############

SUBROUTINE tic(mark)

  INTEGER :: mark
  CALL SYSTEM_CLOCK(mark)
  100 FORMAT('Tic: ', $)
  WRITE(6,100)
  CALL STAMP()

END SUBROUTINE


SUBROUTINE toc(mark)

  INTEGER :: mark
  INTEGER :: now, RATE, CMAX
  REAL :: secs
  INTEGER :: hours, mins

  CALL SYSTEM_CLOCK(now,RATE,CMAX)

  hours = INT((now-mark)/REAL(RATE)/3600.)
  mins = INT((now-mark)/REAL(RATE)/60.) - hours*60
  secs = (now-mark)/REAL(RATE) - hours*3600 - mins*60

  100 FORMAT('Toc: ',$)
  WRITE(*,100)
  CALL STAMP()

  IF ((hours==0).AND.(mins==0)) THEN
    201 FORMAT('Elapsed time: ',F6.3,'s')
    WRITE(*,201) secs
  ELSEIF (hours==0) THEN
    202 FORMAT('Elapsed time: ',I2,'m ',F6.3,'s')
    WRITE(*,202) mins, secs
  ELSE
    203 FORMAT('Elapsed time: ',I4,'h ',I2,'m ',F6.3,'s')
    WRITE(*,203) hours, mins, secs
  END IF

END SUBROUTINE


SUBROUTINE STAMP()

  INTEGER :: now(8)
  CHARACTER(8) :: date
  CHARACTER(10) :: time
  CHARACTER(5) :: zone
  CHARACTER(3) :: month

  CALL DATE_AND_TIME(date,time,zone,now)

  SELECT CASE (date(5:6))
    CASE ('01')
      month = 'Jan'
    CASE ('02')
      month = 'Feb'
    CASE ('03')
      month = 'Mar'
    CASE ('04')
      month = 'Apr'
    CASE ('05')
      month = 'May'
    CASE ('06')
      month = 'Jun'
    CASE ('07')
      month = 'Jul'
    CASE ('08')
      month = 'Aug'
    CASE ('09')
      month = 'Sep'
    CASE ('10')
      month = 'Oct'
    CASE ('11')
      month = 'Nov'
    CASE ('12')
      month = 'Dec'
  END SELECT

  100 FORMAT(A,'/',A,'/',A,' @ ',A,':',A,':',A)
  WRITE(6,100) date(7:8), month, date(1:4), time(1:2), time(3:4), time(5:6)

END SUBROUTINE
