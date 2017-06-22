///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
#ifdef __CINT__


#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;


#pragma link C++ class rpwa::fitResult+;
#pragma link C++ class rpwa::complexMatrix-;


#pragma read                                                                                        \
        sourceClass="TCMatrix"                                                                      \
        source="TMatrixD _re; TMatrixD _im"                                                         \
        version="[1-]"                                                                              \
        targetClass="rpwa::complexMatrix"                                                           \
        target="_size1, _size2, _nmbDataElements, _data"                                            \
        include="TMatrixD.h"                                                                        \
        code="{                                                                                     \
{                                                                                                   \
    _size1 = onfile._re.GetNrows();                                                                 \
    _size2 = onfile._re.GetNcols();                                                                 \
    _nmbDataElements = _size1 * _size2;                                                             \
    _data = new std::complex<double>[_nmbDataElements];                                             \
    unsigned int dataIndex = 0;                                                                     \
    for (unsigned int row = 0; row < _size1; ++row) {                                               \
        for (unsigned int col = 0; col < _size2; ++col) {                                           \
            _data[dataIndex++] = std::complex<double>(onfile._re(row, col), onfile._im(row, col));  \
        }                                                                                           \
    }                                                                                               \
    newObj->readMatrix();                                                                           \
}                                                                                                   \
              }";
#pragma read                                                                                        \
        sourceClass="rpwa::complexMatrix"                                                           \
        source="TMatrixD _re; TMatrixD _im"                                                         \
        version="[1]"                                                                               \
        targetClass="rpwa::complexMatrix"                                                           \
        target="_size1, _size2, _nmbDataElements, _data"                                            \
        include="TMatrixD.h"                                                                        \
        code="{                                                                                     \
{                                                                                                   \
    _size1 = onfile._re.GetNrows();                                                                 \
    _size2 = onfile._re.GetNcols();                                                                 \
    _nmbDataElements = _size1 * _size2;                                                             \
    _data = new std::complex<double>[_nmbDataElements];                                             \
    unsigned int dataIndex = 0;                                                                     \
    for (unsigned int row = 0; row < _size1; ++row)                                                 \
        for (unsigned int col = 0; col < _size2; ++col)                                             \
            _data[dataIndex++] = std::complex<double>(onfile._re(row, col), onfile._im(row, col));  \
    newObj->readMatrix();                                                                           \
}                                                                                                   \
              }";

#pragma read                                                                                        \
        sourceClass="rpwa::fitResult"                                                               \
        source="Double_t _massBinCenter"                                                            \
        version="[-6]"                                                                              \
        targetClass="rpwa::fitResult"                                                               \
        target=""                                                                                   \
        code="{                                                                                     \
{                                                                                                   \
    rpwa::multibinBoundariesType bmr;                                                               \
    bmr[\"mass\"] = rpwa::boundaryType(onfile._massBinCenter, onfile._massBinCenter);               \
    newObj->setMultibinBoundaries(bmr);                                                             \
}                                                                                                   \
              }";

#endif
