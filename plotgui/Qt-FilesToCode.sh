#! /bin/sh

#add names inside the ' ' separated by a space, e.g.: 'name1 name2 name3'
uis='GuiPwaMain'
#mocs='GuiPwaMain GuiWaveTreeModel'

for ui in $uis; do
    uic-qt4 -o ui_$ui.h Qt-DesignerFiles/$ui.ui
done

#for moc in $mocs; do
#    moc-qt4 -o moc_$moc.cxx $moc.h
#done