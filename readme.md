### FinderMaker ###

*Author: Isaac Shivvers*
*Date:   July 2014*

This is an interactive Python code to create nice finder charts
 with offset stars.  To use this, you need a fits image with
 your object clearly detected or the (accurate!) coordinates
 of your object.


Example:

    ipython --pylab
    from findermaker import FinderMaker
    F = FinderMaker( ra='09 55 42.14', dec='+69 40 26.0', name='SN 2014J' )

