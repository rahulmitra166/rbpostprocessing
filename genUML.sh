PYFILES=`find lib -name "*.py"`
echo $PYFILES
for i in $PYFILES
do
    pyreverse -o png $i
    N=`echo $i | sed -e 's/lib\///' -e 's/\.py/.png/'`
    mv classes.png doc/$N
done
