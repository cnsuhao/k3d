// K-3D
// Copyright (c) 2002-2005, Romain Behar
//
// Contact: romainbehar@yahoo.com
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

/** \file
		\author Romain Behar (romainbehar@yahoo.com)
*/

#include <k3dsdk/document_plugin_factory.h>
#include <k3d-i18n-config.h>
#include <k3dsdk/imaterial.h>
#include <k3dsdk/material.h>
#include <k3dsdk/material_sink.h>
#include <k3dsdk/measurement.h>
#include <k3dsdk/legacy_mesh_source.h>
#include <k3dsdk/node.h>
#include <k3dsdk/persistent.h>

namespace libk3dprimitives
{

const double TeapotPoints[306][3] =
{
	{1.4,0.0,2.4},
	{1.4,-0.784,2.4},
	{0.784,-1.4,2.4},
	{0.0,-1.4,2.4},
	{1.3375,0.0,2.53125},
	{1.3375,-0.749,2.53125},
	{0.749,-1.3375,2.53125},
	{0.0,-1.3375,2.53125},
	{1.4375,0.0,2.53125},
	{1.4375,-0.805,2.53125},
	{0.805,-1.4375,2.53125},
	{0.0,-1.4375,2.53125},
	{1.5,0.0,2.4},
	{1.5,-0.84,2.4},
	{0.84,-1.5,2.4},
	{0.0,-1.5,2.4},
	{-0.784,-1.4,2.4},
	{-1.4,-0.784,2.4},
	{-1.4,0.0,2.4},
	{-0.749,-1.3375,2.53125},
	{-1.3375,-0.749,2.53125},
	{-1.3375,0.0,2.53125},
	{-0.805,-1.4375,2.53125},
	{-1.4375,-0.805,2.53125},
	{-1.4375,0.0,2.53125},
	{-0.84,-1.5,2.4},
	{-1.5,-0.84,2.4},
	{-1.5,0.0,2.4},
	{-1.4,0.784,2.4},
	{-0.784,1.4,2.4},
	{0.0,1.4,2.4},
	{-1.3375,0.749,2.53125},
	{-0.749,1.3375,2.53125},
	{0.0,1.3375,2.53125},
	{-1.4375,0.805,2.53125},
	{-0.805,1.4375,2.53125},
	{0.0,1.4375,2.53125},
	{-1.5,0.84,2.4},
	{-0.84,1.5,2.4},
	{0.0,1.5,2.4},
	{0.784,1.4,2.4},
	{1.4,0.784,2.4},
	{0.749,1.3375,2.53125},
	{1.3375,0.749,2.53125},
	{0.805,1.4375,2.53125},
	{1.4375,0.805,2.53125},
	{0.84,1.5,2.4},
	{1.5,0.84,2.4},
	{1.75,0.0,1.875},
	{1.75,-0.98,1.875},
	{0.98,-1.75,1.875},
	{0.0,-1.75,1.875},
	{2.0,0.0,1.35},
	{2.0,-1.12,1.35},
	{1.12,-2.0,1.35},
	{0.0,-2.0,1.35},
	{2.0,0.0,0.9},
	{2.0,-1.12,0.9},
	{1.12,-2.0,0.9},
	{0.0,-2.0,0.9},
	{-0.98,-1.75,1.875},
	{-1.75,-0.98,1.875},
	{-1.75,0.0,1.875},
	{-1.12,-2.0,1.35},
	{-2.0,-1.12,1.35},
	{-2.0,0.0,1.35},
	{-1.12,-2.0,0.9},
	{-2.0,-1.12,0.9},
	{-2.0,0.0,0.9},
	{-1.75,0.98,1.875},
	{-0.98,1.75,1.875},
	{0.0,1.75,1.875},
	{-2.0,1.12,1.35},
	{-1.12,2.0,1.35},
	{0.0,2.0,1.35},
	{-2.0,1.12,0.9},
	{-1.12,2.0,0.9},
	{0.0,2.0,0.9},
	{0.98,1.75,1.875},
	{1.75,0.98,1.875},
	{1.12,2.0,1.35},
	{2.0,1.12,1.35},
	{1.12,2.0,0.9},
	{2.0,1.12,0.9},
	{2.0,0.0,0.45},
	{2.0,-1.12,0.45},
	{1.12,-2.0,0.45},
	{0.0,-2.0,0.45},
	{1.5,0.0,0.225},
	{1.5,-0.84,0.225},
	{0.84,-1.5,0.225},
	{0.0,-1.5,0.225},
	{1.5,0.0,0.15},
	{1.5,-0.84,0.15},
	{0.84,-1.5,0.15},
	{0.0,-1.5,0.15},
	{-1.12,-2.0,0.45},
	{-2.0,-1.12,0.45},
	{-2.0,0.0,0.45},
	{-0.84,-1.5,0.225},
	{-1.5,-0.84,0.225},
	{-1.5,0.0,0.225},
	{-0.84,-1.5,0.15},
	{-1.5,-0.84,0.15},
	{-1.5,0.0,0.15},
	{-2.0,1.12,0.45},
	{-1.12,2.0,0.45},
	{0.0,2.0,0.45},
	{-1.5,0.84,0.225},
	{-0.84,1.5,0.225},
	{0.0,1.5,0.225},
	{-1.5,0.84,0.15},
	{-0.84,1.5,0.15},
	{0.0,1.5,0.15},
	{1.12,2.0,0.45},
	{2.0,1.12,0.45},
	{0.84,1.5,0.225},
	{1.5,0.84,0.225},
	{0.84,1.5,0.15},
	{1.5,0.84,0.15},
	{-1.6,0.0,2.025},
	{-1.6,-0.3,2.025},
	{-1.5,-0.3,2.25},
	{-1.5,0.0,2.25},
	{-2.3,0.0,2.025},
	{-2.3,-0.3,2.025},
	{-2.5,-0.3,2.25},
	{-2.5,0.0,2.25},
	{-2.7,0.0,2.025},
	{-2.7,-0.3,2.025},
	{-3.0,-0.3,2.25},
	{-3.0,0.0,2.25},
	{-2.7,0.0,1.8},
	{-2.7,-0.3,1.8},
	{-3.0,-0.3,1.8},
	{-3.0,0.0,1.8},
	{-1.5,0.3,2.25},
	{-1.6,0.3,2.025},
	{-2.5,0.3,2.25},
	{-2.3,0.3,2.025},
	{-3.0,0.3,2.25},
	{-2.7,0.3,2.025},
	{-3.0,0.3,1.8},
	{-2.7,0.3,1.8},
	{-2.7,0.0,1.575},
	{-2.7,-0.3,1.575},
	{-3.0,-0.3,1.35},
	{-3.0,0.0,1.35},
	{-2.5,0.0,1.125},
	{-2.5,-0.3,1.125},
	{-2.65,-0.3,0.9375},
	{-2.65,0.0,0.9375},
	{-2.0,-0.3,0.9},
	{-1.9,-0.3,0.6},
	{-1.9,0.0,0.6},
	{-3.0,0.3,1.35},
	{-2.7,0.3,1.575},
	{-2.65,0.3,0.9375},
	{-2.5,0.3,1.125},
	{-1.9,0.3,0.6},
	{-2.0,0.3,0.9},
	{1.7,0.0,1.425},
	{1.7,-0.66,1.425},
	{1.7,-0.66,0.6},
	{1.7,0.0,0.6},
	{2.6,0.0,1.425},
	{2.6,-0.66,1.425},
	{3.1,-0.66,0.825},
	{3.1,0.0,0.825},
	{2.3,0.0,2.1},
	{2.3,-0.25,2.1},
	{2.4,-0.25,2.025},
	{2.4,0.0,2.025},
	{2.7,0.0,2.4},
	{2.7,-0.25,2.4},
	{3.3,-0.25,2.4},
	{3.3,0.0,2.4},
	{1.7,0.66,0.6},
	{1.7,0.66,1.425},
	{3.1,0.66,0.825},
	{2.6,0.66,1.425},
	{2.4,0.25,2.025},
	{2.3,0.25,2.1},
	{3.3,0.25,2.4},
	{2.7,0.25,2.4},
	{2.8,0.0,2.475},
	{2.8,-0.25,2.475},
	{3.525,-0.25,2.49375},
	{3.525,0.0,2.49375},
	{2.9,0.0,2.475},
	{2.9,-0.15,2.475},
	{3.45,-0.15,2.5125},
	{3.45,0.0,2.5125},
	{2.8,0.0,2.4},
	{2.8,-0.15,2.4},
	{3.2,-0.15,2.4},
	{3.2,0.0,2.4},
	{3.525,0.25,2.49375},
	{2.8,0.25,2.475},
	{3.45,0.15,2.5125},
	{2.9,0.15,2.475},
	{3.2,0.15,2.4},
	{2.8,0.15,2.4},
	{0.0,0.0,3.15},
	{0.0,-0.002,3.15},
	{0.002,0.0,3.15},
	{0.8,0.0,3.15},
	{0.8,-0.45,3.15},
	{0.45,-0.8,3.15},
	{0.0,-0.8,3.15},
	{0.0,0.0,2.85},
	{0.2,0.0,2.7},
	{0.2,-0.112,2.7},
	{0.112,-0.2,2.7},
	{0.0,-0.2,2.7},
	{-0.002,0.0,3.15},
	{-0.45,-0.8,3.15},
	{-0.8,-0.45,3.15},
	{-0.8,0.0,3.15},
	{-0.112,-0.2,2.7},
	{-0.2,-0.112,2.7},
	{-0.2,0.0,2.7},
	{0.0,0.002,3.15},
	{-0.8,0.45,3.15},
	{-0.45,0.8,3.15},
	{0.0,0.8,3.15},
	{-0.2,0.112,2.7},
	{-0.112,0.2,2.7},
	{0.0,0.2,2.7},
	{0.45,0.8,3.15},
	{0.8,0.45,3.15},
	{0.112,0.2,2.7},
	{0.2,0.112,2.7},
	{0.4,0.0,2.55},
	{0.4,-0.224,2.55},
	{0.224,-0.4,2.55},
	{0.0,-0.4,2.55},
	{1.3,0.0,2.55},
	{1.3,-0.728,2.55},
	{0.728,-1.3,2.55},
	{0.0,-1.3,2.55},
	{1.3,0.0,2.4},
	{1.3,-0.728,2.4},
	{0.728,-1.3,2.4},
	{0.0,-1.3,2.4},
	{-0.224,-0.4,2.55},
	{-0.4,-0.224,2.55},
	{-0.4,0.0,2.55},
	{-0.728,-1.3,2.55},
	{-1.3,-0.728,2.55},
	{-1.3,0.0,2.55},
	{-0.728,-1.3,2.4},
	{-1.3,-0.728,2.4},
	{-1.3,0.0,2.4},
	{-0.4,0.224,2.55},
	{-0.224,0.4,2.55},
	{0.0,0.4,2.55},
	{-1.3,0.728,2.55},
	{-0.728,1.3,2.55},
	{0.0,1.3,2.55},
	{-1.3,0.728,2.4},
	{-0.728,1.3,2.4},
	{0.0,1.3,2.4},
	{0.224,0.4,2.55},
	{0.4,0.224,2.55},
	{0.728,1.3,2.55},
	{1.3,0.728,2.55},
	{0.728,1.3,2.4},
	{1.3,0.728,2.4},
	{0.0,0.0,0.0},
	{1.5,0.0,0.15},
	{1.5,0.84,0.15},
	{0.84,1.5,0.15},
	{0.0,1.5,0.15},
	{1.5,0.0,0.075},
	{1.5,0.84,0.075},
	{0.84,1.5,0.075},
	{0.0,1.5,0.075},
	{1.425,0.0,0.0},
	{1.425,0.798,0.0},
	{0.798,1.425,0.0},
	{0.0,1.425,0.0},
	{-0.84,1.5,0.15},
	{-1.5,0.84,0.15},
	{-1.5,0.0,0.15},
	{-0.84,1.5,0.075},
	{-1.5,0.84,0.075},
	{-1.5,0.0,0.075},
	{-0.798,1.425,0.0},
	{-1.425,0.798,0.0},
	{-1.425,0.0,0.0},
	{-1.5,-0.84,0.15},
	{-0.84,-1.5,0.15},
	{0.0,-1.5,0.15},
	{-1.5,-0.84,0.075},
	{-0.84,-1.5,0.075},
	{0.0,-1.5,0.075},
	{-1.425,-0.798,0.0},
	{-0.798,-1.425,0.0},
	{0.0,-1.425,0.0},
	{0.84,-1.5,0.15},
	{1.5,-0.84,0.15},
	{0.84,-1.5,0.075},
	{1.5,-0.84,0.075},
	{0.798,-1.425,0.0},
	{1.425,-0.798,0.0}
};

const unsigned long TeapotPatches[32][16] =
{
	// Rim
	{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16},
	{4,17,18,19,8,20,21,22,12,23,24,25,16,26,27,28},
	{19,29,30,31,22,32,33,34,25,35,36,37,28,38,39,40},
	{31,41,42,1,34,43,44,5,37,45,46,9,40,47,48,13},
	// Body
	{13,14,15,16,49,50,51,52,53,54,55,56,57,58,59,60},
	{16,26,27,28,52,61,62,63,56,64,65,66,60,67,68,69},
	{28,38,39,40,63,70,71,72,66,73,74,75,69,76,77,78},
	{40,47,48,13,72,79,80,49,75,81,82,53,78,83,84,57},
	{57,58,59,60,85,86,87,88,89,90,91,92,93,94,95,96},
	{60,67,68,69,88,97,98,99,92,100,101,102,96,103,104,105},
	{69,76,77,78,99,106,107,108,102,109,110,111,105,112,113,114},
	{78,83,84,57,108,115,116,85,111,117,118,89,114,119,120,93},
	// Handle
	{121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136},
	{124,137,138,121,128,139,140,125,132,141,142,129,136,143,144,133},
	{133,134,135,136,145,146,147,148,149,150,151,152,69,153,154,155},
	{136,143,144,133,148,156,157,145,152,158,159,149,155,160,161,69},
	// Spout
	{162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177},
	{165,178,179,162,169,180,181,166,173,182,183,170,177,184,185,174},
	{174,175,176,177,186,187,188,189,190,191,192,193,194,195,196,197},
	{177,184,185,174,189,198,199,186,193,200,201,190,197,202,203,194},
	// Lid
	{204,204,204,204,207,208,209,210,211,211,211,211,212,213,214,215},
	{204,204,204,204,210,217,218,219,211,211,211,211,215,220,221,222},
	{204,204,204,204,219,224,225,226,211,211,211,211,222,227,228,229},
	{204,204,204,204,226,230,231,207,211,211,211,211,229,232,233,212},
	{212,213,214,215,234,235,236,237,238,239,240,241,242,243,244,245},
	{215,220,221,222,237,246,247,248,241,249,250,251,245,252,253,254},
	{222,227,228,229,248,255,256,257,251,258,259,260,254,261,262,263},
	{229,232,233,212,257,264,265,234,260,266,267,238,263,268,269,242},
	// Bottom (not on original teapot)
	{270,270,270,270,279,280,281,282,275,276,277,278,271,272,273,274},
	{270,270,270,270,282,289,290,291,278,286,287,288,274,283,284,285},
	{270,270,270,270,291,298,299,300,288,295,296,297,285,292,293,294},
	{270,270,270,270,300,305,306,279,297,303,304,275,294,301,302,271}
};

const double TeacupPoints[251][3] =
{
	{0.409091,0.772727,0.0},
	{0.409091,0.772727,-0.229091},
	{0.229091,0.772727,-0.409091},
	{0.0,0.772727,-0.409091},
	{0.409091,0.886364,0.0},
	{0.409091,0.886364,-0.229091},
	{0.229091,0.886364,-0.409091},
	{0.0,0.886364,-0.409091},
	{0.454545,0.886364,0.0},
	{0.454545,0.886364,-0.254545},
	{0.254545,0.886364,-0.454545},
	{0.0,0.886364,-0.454545},
	{0.454545,0.772727,0.0},
	{0.454545,0.772727,-0.254545},
	{0.254545,0.772727,-0.454545},
	{0.0,0.772727,-0.454545},
	{-0.229091,0.772727,-0.409091},
	{-0.409091,0.772727,-0.229091},
	{-0.409091,0.772727,0.0},
	{-0.229091,0.886364,-0.409091},
	{-0.409091,0.886364,-0.229091},
	{-0.409091,0.886364,0.0},
	{-0.254545,0.886364,-0.454545},
	{-0.454545,0.886364,-0.254545},
	{-0.454545,0.886364,0.0},
	{-0.254545,0.772727,-0.454545},
	{-0.454545,0.772727,-0.254545},
	{-0.454545,0.772727,0.0},
	{-0.409091,0.772727,0.229091},
	{-0.229091,0.772727,0.409091},
	{0.0,0.772727,0.409091},
	{-0.409091,0.886364,0.229091},
	{-0.229091,0.886364,0.409091},
	{0.0,0.886364,0.409091},
	{-0.454545,0.886364,0.254545},
	{-0.254545,0.886364,0.454545},
	{0.0,0.886364,0.454545},
	{-0.454545,0.772727,0.254545},
	{-0.254545,0.772727,0.454545},
	{0.0,0.772727,0.454545},
	{0.229091,0.772727,0.409091},
	{0.409091,0.772727,0.229091},
	{0.229091,0.886364,0.409091},
	{0.409091,0.886364,0.229091},
	{0.254545,0.886364,0.454545},
	{0.454545,0.886364,0.254545},
	{0.254545,0.772727,0.454545},
	{0.454545,0.772727,0.254545},
	{0.454545,0.545455,0.0},
	{0.454545,0.545455,-0.254545},
	{0.254545,0.545455,-0.454545},
	{0.0,0.545455,-0.454545},
	{0.454545,0.272727,0.0},
	{0.454545,0.272727,-0.254545},
	{0.254545,0.272727,-0.454545},
	{0.0,0.272727,-0.454545},
	{0.318182,0.0454545,0.0},
	{0.318182,0.0454545,-0.178182},
	{0.178182,0.0454545,-0.318182},
	{0.0,0.0454545,-0.318182},
	{-0.254545,0.545455,-0.454545},
	{-0.454545,0.545455,-0.254545},
	{-0.454545,0.545455,0.0},
	{-0.254545,0.272727,-0.454545},
	{-0.454545,0.272727,-0.254545},
	{-0.454545,0.272727,0.0},
	{-0.178182,0.0454545,-0.318182},
	{-0.318182,0.0454545,-0.178182},
	{-0.318182,0.0454545,0.0},
	{-0.454545,0.545455,0.254545},
	{-0.254545,0.545455,0.454545},
	{0.0,0.545455,0.454545},
	{-0.454545,0.272727,0.254545},
	{-0.254545,0.272727,0.454545},
	{0.0,0.272727,0.454545},
	{-0.318182,0.0454545,0.178182},
	{-0.178182,0.0454545,0.318182},
	{0.0,0.0454545,0.318182},
	{0.254545,0.545455,0.454545},
	{0.454545,0.545455,0.254545},
	{0.254545,0.272727,0.454545},
	{0.454545,0.272727,0.254545},
	{0.178182,0.0454545,0.318182},
	{0.318182,0.0454545,0.178182},
	{0.545455,0.0454545,0.0},
	{0.545455,0.0454545,-0.305455},
	{0.305455,0.0454545,-0.545455},
	{0.0,0.0454545,-0.545455},
	{0.727273,0.136364,0.0},
	{0.727273,0.136364,-0.407273},
	{0.407273,0.136364,-0.727273},
	{0.0,0.136364,-0.727273},
	{0.909091,0.136364,0.0},
	{0.909091,0.136364,-0.509091},
	{0.509091,0.136364,-0.909091},
	{0.0,0.136364,-0.909091},
	{-0.305455,0.0454545,-0.545455},
	{-0.545455,0.0454545,-0.305455},
	{-0.545455,0.0454545,0.0},
	{-0.407273,0.136364,-0.727273},
	{-0.727273,0.136364,-0.407273},
	{-0.727273,0.136364,0.0},
	{-0.509091,0.136364,-0.909091},
	{-0.909091,0.136364,-0.509091},
	{-0.909091,0.136364,0.0},
	{-0.545455,0.0454545,0.305455},
	{-0.305455,0.0454545,0.545455},
	{0.0,0.0454545,0.545455},
	{-0.727273,0.136364,0.407273},
	{-0.407273,0.136364,0.727273},
	{0.0,0.136364,0.727273},
	{-0.909091,0.136364,0.509091},
	{-0.509091,0.136364,0.909091},
	{0.0,0.136364,0.909091},
	{0.305455,0.0454545,0.545455},
	{0.545455,0.0454545,0.305455},
	{0.407273,0.136364,0.727273},
	{0.727273,0.136364,0.407273},
	{0.509091,0.136364,0.909091},
	{0.909091,0.136364,0.509091},
	{1.0,0.136364,0.0},
	{1.0,0.136364,-0.56},
	{0.56,0.136364,-1.0},
	{0.0,0.136364,-1.0},
	{1.0,0.0909091,0.0},
	{1.0,0.0909091,-0.56},
	{0.56,0.0909091,-1.0},
	{0.0,0.0909091,-1.0},
	{0.909091,0.0909091,0.0},
	{0.909091,0.0909091,-0.509091},
	{0.509091,0.0909091,-0.909091},
	{0.0,0.0909091,-0.909091},
	{-0.56,0.136364,-1.0},
	{-1.0,0.136364,-0.56},
	{-1.0,0.136364,0.0},
	{-0.56,0.0909091,-1.0},
	{-1.0,0.0909091,-0.56},
	{-1.0,0.0909091,0.0},
	{-0.509091,0.0909091,-0.909091},
	{-0.909091,0.0909091,-0.509091},
	{-0.909091,0.0909091,0.0},
	{-1.0,0.136364,0.56},
	{-0.56,0.136364,1.0},
	{0.0,0.136364,1.0},
	{-1.0,0.0909091,0.56},
	{-0.56,0.0909091,1.0},
	{0.0,0.0909091,1.0},
	{-0.909091,0.0909091,0.509091},
	{-0.509091,0.0909091,0.909091},
	{0.0,0.0909091,0.909091},
	{0.56,0.136364,1.0},
	{1.0,0.136364,0.56},
	{0.56,0.0909091,1.0},
	{1.0,0.0909091,0.56},
	{0.509091,0.0909091,0.909091},
	{0.909091,0.0909091,0.509091},
	{0.727273,0.0909091,0.0},
	{0.727273,0.0909091,-0.407273},
	{0.407273,0.0909091,-0.727273},
	{0.0,0.0909091,-0.727273},
	{0.545455,0.0,0.0},
	{0.545455,0.0,-0.305455},
	{0.305455,0.0,-0.545455},
	{0.0,0.0,-0.545455},
	{0.318182,0.0,0.0},
	{0.318182,0.0,-0.178182},
	{0.178182,0.0,-0.318182},
	{0.0,0.0,-0.318182},
	{-0.407273,0.0909091,-0.727273},
	{-0.727273,0.0909091,-0.407273},
	{-0.727273,0.0909091,0.0},
	{-0.305455,0.0,-0.545455},
	{-0.545455,0.0,-0.305455},
	{-0.545455,0.0,0.0},
	{-0.178182,0.0,-0.318182},
	{-0.318182,0.0,-0.178182},
	{-0.318182,0.0,0.0},
	{-0.727273,0.0909091,0.407273},
	{-0.407273,0.0909091,0.727273},
	{0.0,0.0909091,0.727273},
	{-0.545455,0.0,0.305455},
	{-0.305455,0.0,0.545455},
	{0.0,0.0,0.545455},
	{-0.318182,0.0,0.178182},
	{-0.178182,0.0,0.318182},
	{0.0,0.0,0.318182},
	{0.407273,0.0909091,0.727273},
	{0.727273,0.0909091,0.407273},
	{0.305455,0.0,0.545455},
	{0.545455,0.0,0.305455},
	{0.178182,0.0,0.318182},
	{0.318182,0.0,0.178182},
	{0.272727,0.0454545,0.0},
	{0.272727,0.0454545,-0.152727},
	{0.152727,0.0454545,-0.272727},
	{0.0,0.0454545,-0.272727},
	{0.409091,0.272727,0.0},
	{0.409091,0.272727,-0.229091},
	{0.229091,0.272727,-0.409091},
	{0.0,0.272727,-0.409091},
	{0.409091,0.545455,0.0},
	{0.409091,0.545455,-0.229091},
	{0.229091,0.545455,-0.409091},
	{0.0,0.545455,-0.409091},
	{-0.152727,0.0454545,-0.272727},
	{-0.272727,0.0454545,-0.152727},
	{-0.272727,0.0454545,0.0},
	{-0.229091,0.272727,-0.409091},
	{-0.409091,0.272727,-0.229091},
	{-0.409091,0.272727,0.0},
	{-0.229091,0.545455,-0.409091},
	{-0.409091,0.545455,-0.229091},
	{-0.409091,0.545455,0.0},
	{-0.272727,0.0454545,0.152727},
	{-0.152727,0.0454545,0.272727},
	{0.0,0.0454545,0.272727},
	{-0.409091,0.272727,0.229091},
	{-0.229091,0.272727,0.409091},
	{0.0,0.272727,0.409091},
	{-0.409091,0.545455,0.229091},
	{-0.229091,0.545455,0.409091},
	{0.0,0.545455,0.409091},
	{0.152727,0.0454545,0.272727},
	{0.272727,0.0454545,0.152727},
	{0.229091,0.272727,0.409091},
	{0.409091,0.272727,0.229091},
	{0.229091,0.545455,0.409091},
	{0.409091,0.545455,0.229091},
	{-0.454545,0.704545,0.0},
	{-0.454545,0.704545,-0.0454545},
	{-0.454545,0.772727,-0.0454545},
	{-0.772727,0.863636,0.0},
	{-0.772727,0.863636,-0.0454545},
	{-0.818182,0.954545,-0.0454545},
	{-0.818182,0.954545,0.0},
	{-0.772727,0.522727,0.0},
	{-0.772727,0.522727,-0.0454545},
	{-0.909091,0.477273,-0.0454545},
	{-0.909091,0.477273,0.0},
	{-0.409091,0.363636,0.0},
	{-0.409091,0.363636,-0.0454545},
	{-0.409091,0.295455,-0.0454545},
	{-0.409091,0.295455,0.0},
	{-0.454545,0.772727,0.0454545},
	{-0.454545,0.704545,0.0454545},
	{-0.818182,0.954545,0.0454545},
	{-0.772727,0.863636,0.0454545},
	{-0.909091,0.477273,0.0454545},
	{-0.772727,0.522727,0.0454545},
	{-0.409091,0.295455,0.0454545},
	{-0.409091,0.363636,0.0454545}
};

const unsigned long TeacupPatches[26][16] =
{
	{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16},
	{4,17,18,19,8,20,21,22,12,23,24,25,16,26,27,28},
	{19,29,30,31,22,32,33,34,25,35,36,37,28,38,39,40},
	{31,41,42,1,34,43,44,5,37,45,46,9,40,47,48,13},
	{13,14,15,16,49,50,51,52,53,54,55,56,57,58,59,60},
	{16,26,27,28,52,61,62,63,56,64,65,66,60,67,68,69},
	{28,38,39,40,63,70,71,72,66,73,74,75,69,76,77,78},
	{40,47,48,13,72,79,80,49,75,81,82,53,78,83,84,57},
	{193,194,195,196,197,198,199,200,201,202,203,204,1,2,3,4},
	{196,205,206,207,200,208,209,210,204,211,212,213,4,17,18,19},
	{207,214,215,216,210,217,218,219,213,220,221,222,19,29,30,31},
	{216,223,224,193,219,225,226,197,222,227,228,201,31,41,42,1},
	{229,230,231,28,232,233,234,235,236,237,238,239,240,241,242,243},
	{28,244,245,229,235,246,247,232,239,248,249,236,243,250,251,240},
	{57,58,59,60,85,86,87,88,89,90,91,92,93,94,95,96},
	{60,67,68,69,88,97,98,99,92,100,101,102,96,103,104,105},
	{69,76,77,78,99,106,107,108,102,109,110,111,105,112,113,114},
	{78,83,84,57,108,115,116,85,111,117,118,89,114,119,120,93},
	{93,94,95,96,121,122,123,124,125,126,127,128,129,130,131,132},
	{96,103,104,105,124,133,134,135,128,136,137,138,132,139,140,141},
	{105,112,113,114,135,142,143,144,138,145,146,147,141,148,149,150},
	{114,119,120,93,144,151,152,121,147,153,154,125,150,155,156,129},
	{129,130,131,132,157,158,159,160,161,162,163,164,165,166,167,168},
	{132,139,140,141,160,169,170,171,164,172,173,174,168,175,176,177},
	{141,148,149,150,171,178,179,180,174,181,182,183,177,184,185,186},
	{150,155,156,129,180,187,188,157,183,189,190,161,186,191,192,165}
};

const double TeaspoonPoints[256][3] =
{
	{-0.000107143,0.205357,0.0},
	{0.0,0.196429,-0.0178571},
	{0.0,0.196429,-0.0178571},
	{0.000107143,0.205357,0.0},
	{-0.0535714,0.205357,0.0},
	{-0.0222714,0.178571,-0.0534286},
	{0.0222714,0.178571,-0.0534286},
	{0.0535714,0.205357,0.0},
	{-0.107143,0.0952429,-0.0178571},
	{-0.0446429,0.0952429,-0.0892857},
	{0.0446429,0.0952429,-0.0892857},
	{0.107143,0.0952429,-0.0178571},
	{-0.107143,0.0,-0.0178571},
	{-0.0446429,0.0,-0.0892857},
	{0.0446429,0.0,-0.0892857},
	{0.107143,0.0,-0.0178571},
	{0.000107143,0.205357,0.0},
	{0.000135714,0.207589,0.00446429},
	{0.000157143,0.216518,0.00446429},
	{0.000125,0.214286,0.0},
	{0.0535714,0.205357,0.0},
	{0.0613964,0.212054,0.0133571},
	{0.0714286,0.220982,0.015625},
	{0.0625,0.214286,0.0},
	{0.107143,0.0952429,-0.0178571},
	{0.122768,0.0952429,0.0},
	{0.142857,0.0952429,0.00446429},
	{0.125,0.0952429,-0.0178571},
	{0.107143,0.0,-0.0178571},
	{0.122768,0.0,0.0},
	{0.142857,0.0,0.00446429},
	{0.125,0.0,-0.0178571},
	{0.000125,0.214286,0.0},
	{0.0,0.205357,-0.0178571},
	{0.0,0.205357,-0.0178571},
	{-0.000125,0.214286,0.0},
	{0.0625,0.214286,0.0},
	{0.0267857,0.1875,-0.0625},
	{-0.0267857,0.1875,-0.0625},
	{-0.0625,0.214286,0.0},
	{0.125,0.0952429,-0.0178571},
	{0.0535714,0.0952429,-0.107143},
	{-0.0535714,0.0952429,-0.107143},
	{-0.125,0.0952429,-0.0178571},
	{0.125,0.0,-0.0178571},
	{0.0535714,0.0,-0.107143},
	{-0.0535714,0.0,-0.107143},
	{-0.125,0.0,-0.0178571},
	{-0.000125,0.214286,0.0},
	{-0.000157143,0.216518,0.00446429},
	{-0.000135714,0.207589,0.00446429},
	{-0.000107143,0.205357,0.0},
	{-0.0625,0.214286,0.0},
	{-0.0714286,0.220982,0.015625},
	{-0.0613964,0.212054,0.0133571},
	{-0.0535714,0.205357,0.0},
	{-0.125,0.0952429,-0.0178571},
	{-0.142857,0.0952429,0.00446429},
	{-0.122768,0.0952429,0.0},
	{-0.107143,0.0952429,-0.0178571},
	{-0.125,0.0,-0.0178571},
	{-0.142857,0.0,0.00446429},
	{-0.122768,0.0,0.0},
	{-0.107143,0.0,-0.0178571},
	{-0.107143,0.0,-0.0178571},
	{-0.0446429,0.0,-0.0892857},
	{0.0446429,0.0,-0.0892857},
	{0.107143,0.0,-0.0178571},
	{-0.107143,-0.142857,-0.0178571},
	{-0.0446429,-0.142857,-0.0892857},
	{0.0446429,-0.142857,-0.0892857},
	{0.107143,-0.142857,-0.0178571},
	{-0.0133929,-0.160714,0.0386893},
	{-0.00557857,-0.160714,0.0386893},
	{0.00557857,-0.160714,0.0386893},
	{0.0133929,-0.160714,0.0386893},
	{-0.0133929,-0.25,0.0535714},
	{-0.00557857,-0.25,0.0535714},
	{0.00557857,-0.25,0.0535714},
	{0.0133929,-0.25,0.0535714},
	{0.107143,0.0,-0.0178571},
	{0.122768,0.0,0.0},
	{0.142857,0.0,0.00446429},
	{0.125,0.0,-0.0178571},
	{0.107143,-0.142857,-0.0178571},
	{0.122768,-0.142857,0.0},
	{0.142857,-0.142857,0.00446429},
	{0.125,-0.142857,-0.0178571},
	{0.0133929,-0.160714,0.0386893},
	{0.0153464,-0.160714,0.0386893},
	{0.0178571,-0.160714,0.0314357},
	{0.015625,-0.160714,0.0297607},
	{0.0133929,-0.25,0.0535714},
	{0.0153464,-0.25,0.0535714},
	{0.0178571,-0.25,0.0463179},
	{0.015625,-0.25,0.0446429},
	{0.125,0.0,-0.0178571},
	{0.0535714,0.0,-0.107143},
	{-0.0535714,0.0,-0.107143},
	{-0.125,0.0,-0.0178571},
	{0.125,-0.142857,-0.0178571},
	{0.0535714,-0.142857,-0.107143},
	{-0.0535714,-0.142857,-0.107143},
	{-0.125,-0.142857,-0.0178571},
	{0.015625,-0.160714,0.0297607},
	{0.00669643,-0.160714,0.0230643},
	{-0.00781071,-0.160714,0.0208321},
	{-0.015625,-0.160714,0.0297607},
	{0.015625,-0.25,0.0446429},
	{0.00669643,-0.25,0.0379464},
	{-0.00781071,-0.25,0.0357143},
	{-0.015625,-0.25,0.0446429},
	{-0.125,0.0,-0.0178571},
	{-0.142857,0.0,0.00446429},
	{-0.122768,0.0,0.0},
	{-0.107143,0.0,-0.0178571},
	{-0.125,-0.142857,-0.0178571},
	{-0.142857,-0.142857,0.00446429},
	{-0.122768,-0.142857,0.0},
	{-0.107143,-0.142857,-0.0178571},
	{-0.015625,-0.160714,0.0297607},
	{-0.0175786,-0.160714,0.0319929},
	{-0.0153464,-0.160714,0.0386893},
	{-0.0133929,-0.160714,0.0386893},
	{-0.015625,-0.25,0.0446429},
	{-0.0175786,-0.25,0.046875},
	{-0.0153464,-0.25,0.0535714},
	{-0.0133929,-0.25,0.0535714},
	{-0.0133929,-0.25,0.0535714},
	{-0.00557857,-0.25,0.0535714},
	{0.00557857,-0.25,0.0535714},
	{0.0133929,-0.25,0.0535714},
	{-0.0133929,-0.46425,0.0892857},
	{-0.00557857,-0.46425,0.0892857},
	{0.00557857,-0.46425,0.0892857},
	{0.0133929,-0.46425,0.0892857},
	{-0.0446429,-0.678571,0.0535714},
	{-0.00892857,-0.678571,0.0625},
	{0.00892857,-0.678571,0.0625},
	{0.0446429,-0.678571,0.0535714},
	{-0.0446429,-0.857143,0.0357143},
	{-0.00892857,-0.857143,0.0446429},
	{0.00892857,-0.857143,0.0446429},
	{0.0446429,-0.857143,0.0357143},
	{0.0133929,-0.25,0.0535714},
	{0.0153464,-0.25,0.0535714},
	{0.0178571,-0.25,0.0463179},
	{0.015625,-0.25,0.0446429},
	{0.0133929,-0.46425,0.0892857},
	{0.0153464,-0.464286,0.0892857},
	{0.0178571,-0.46425,0.0820321},
	{0.015625,-0.46425,0.0803571},
	{0.0446429,-0.678571,0.0535714},
	{0.0535714,-0.678571,0.0513393},
	{0.0535714,-0.678571,0.0334821},
	{0.0446429,-0.678571,0.0357143},
	{0.0446429,-0.857143,0.0357143},
	{0.0535714,-0.857143,0.0334821},
	{0.0535714,-0.857143,0.015625},
	{0.0446429,-0.857143,0.0178571},
	{0.015625,-0.25,0.0446429},
	{0.00669643,-0.25,0.0379464},
	{-0.00781071,-0.25,0.0357143},
	{-0.015625,-0.25,0.0446429},
	{0.015625,-0.46425,0.0803571},
	{0.00669643,-0.464286,0.0736607},
	{-0.00781071,-0.46425,0.0714286},
	{-0.015625,-0.46425,0.0803571},
	{0.0446429,-0.678571,0.0357143},
	{0.00892857,-0.678571,0.0446429},
	{-0.00892857,-0.678571,0.0446429},
	{-0.0446429,-0.678571,0.0357143},
	{0.0446429,-0.857143,0.0178571},
	{0.00892857,-0.857143,0.0267857},
	{-0.00892857,-0.857143,0.0267857},
	{-0.0446429,-0.857143,0.0178571},
	{-0.015625,-0.25,0.0446429},
	{-0.0175786,-0.25,0.046875},
	{-0.0153464,-0.25,0.0535714},
	{-0.0133929,-0.25,0.0535714},
	{-0.015625,-0.46425,0.0803571},
	{-0.0175786,-0.464286,0.0825893},
	{-0.0153464,-0.464286,0.0892857},
	{-0.0133929,-0.46425,0.0892857},
	{-0.0446429,-0.678571,0.0357143},
	{-0.0535714,-0.678571,0.0334821},
	{-0.0535714,-0.678571,0.0513393},
	{-0.0446429,-0.678571,0.0535714},
	{-0.0446429,-0.857143,0.0178571},
	{-0.0535714,-0.857143,0.015625},
	{-0.0535714,-0.857143,0.0334821},
	{-0.0446429,-0.857143,0.0357143},
	{-0.0446429,-0.857143,0.0357143},
	{-0.00892857,-0.857143,0.0446429},
	{0.00892857,-0.857143,0.0446429},
	{0.0446429,-0.857143,0.0357143},
	{-0.0446429,-0.928571,0.0285714},
	{-0.00892857,-0.928571,0.0375},
	{0.00892857,-0.928571,0.0375},
	{0.0446429,-0.928571,0.0285714},
	{-0.0539286,-0.999643,0.0178571},
	{0.000357143,-0.999643,0.0178571},
	{0.0,-0.999643,0.0178571},
	{0.0535714,-0.999643,0.0178571},
	{-0.000357143,-1,0.0178571},
	{0.000357143,-1,0.0178571},
	{0.0,-1,0.0178571},
	{0.0,-1,0.0178571},
	{0.0446429,-0.857143,0.0357143},
	{0.0535714,-0.857143,0.0334821},
	{0.0535714,-0.857143,0.015625},
	{0.0446429,-0.857143,0.0178571},
	{0.0446429,-0.928571,0.0285714},
	{0.0535714,-0.928571,0.0263393},
	{0.0535714,-0.928571,0.00848214},
	{0.0446429,-0.928571,0.0107143},
	{0.0535714,-0.999643,0.0178571},
	{0.0669643,-0.999643,0.0178571},
	{0.0673214,-0.999643,0.0},
	{0.0539286,-0.999643,0.0},
	{0.0,-1,0.0178571},
	{0.0,-1,0.0178571},
	{0.000357143,-1,0.0},
	{0.000357143,-1,0.0},
	{0.0446429,-0.857143,0.0178571},
	{0.00892857,-0.857143,0.0267857},
	{-0.00892857,-0.857143,0.0267857},
	{-0.0446429,-0.857143,0.0178571},
	{0.0446429,-0.928571,0.0107143},
	{0.00892857,-0.928571,0.0196429},
	{-0.00892857,-0.928571,0.0196429},
	{-0.0446429,-0.928571,0.0107143},
	{0.0539286,-0.999643,0.0},
	{0.000357143,-0.999643,0.0},
	{-0.000357143,-0.999643,0.0},
	{-0.0539286,-0.999643,0.0},
	{0.000357143,-1,0.0},
	{0.000357143,-1,0.0},
	{-0.000357143,-1,0.0},
	{-0.000357143,-1,0.0},
	{-0.0446429,-0.857143,0.0178571},
	{-0.0535714,-0.857143,0.015625},
	{-0.0535714,-0.857143,0.0334821},
	{-0.0446429,-0.857143,0.0357143},
	{-0.0446429,-0.928571,0.0107143},
	{-0.0535714,-0.928571,0.00848214},
	{-0.0535714,-0.928571,0.0263393},
	{-0.0446429,-0.928571,0.0285714},
	{-0.0539286,-0.999643,0.0},
	{-0.0673214,-0.999643,0.0},
	{-0.0675,-0.999643,0.0178571},
	{-0.0539286,-0.999643,0.0178571},
	{-0.000357143,-1,0.0},
	{-0.000357143,-1,0.0},
	{-0.000535714,-1,0.0178571},
	{-0.000357143,-1,0.0178571}
};

const unsigned long TeaspoonPatches[16][16] =
{
	{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16},
	{17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32},
	{33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48},
	{49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64},
	{65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80},
	{81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96},
	{97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112},
	{113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128},
	{129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144},
	{145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160},
	{161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176},
	{177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192},
	{193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208},
	{209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224},
	{225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240},
	{241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256}
};

/////////////////////////////////////////////////////////////////////////////
// newell_primitive_implementation

class newell_primitive_implementation :
	public k3d::material_sink<k3d::legacy::mesh_source<k3d::persistent<k3d::node> > >
{
	typedef k3d::material_sink<k3d::legacy::mesh_source<k3d::persistent<k3d::node> > > base;

public:
	newell_primitive_implementation(k3d::iplugin_factory& Factory, k3d::idocument& Document) :
		base(Factory, Document),
		m_type(init_owner(*this) + init_name("type") + init_label(_("Primitive")) + init_description(_("Primitive type (teapot, teacup or teaspoon")) + init_value(TEAPOT) + init_enumeration(type_values())),
		m_size(init_owner(*this) + init_name("size") + init_label(_("Size")) + init_description(_("Primitive size (scale)")) + init_value(1.0) + init_step_increment(0.1) + init_units(typeid(k3d::measurement::scalar)))
	{
		m_material.changed_signal().connect(make_reset_mesh_slot());
		m_type.changed_signal().connect(make_reset_mesh_slot());
		m_size.changed_signal().connect(make_reset_mesh_slot());
	}

	void on_initialize_mesh(k3d::legacy::mesh& Mesh)
	{
		const double size = m_size.pipeline_value();
		k3d::imaterial* const material = m_material.pipeline_value();

		// Creates Newell's primitives using his bezier patches ...
		switch(m_type.pipeline_value())
		{
			case TEAPOT:
			{
				// Min : -3.0, max : 3.525 -> 20/6.525 = 3.065
				for(unsigned long i = 0; i < 306; i++)
					Mesh.points.push_back(new k3d::legacy::point(size * k3d::point3(TeapotPoints[i][0], TeapotPoints[i][1], TeapotPoints[i][2])));

				// Rim
				for(unsigned long i = 0; i < 4; i++)
				{
					k3d::legacy::bicubic_patch* const patch = new k3d::legacy::bicubic_patch();
					patch->material = material;

					for(unsigned long j = 0; j < 16; j++)
						patch->control_points[j] = Mesh.points[TeapotPatches[i][j]-1];

					Mesh.bicubic_patches.push_back(patch);
				}
				// Body
				for(unsigned long i = 4; i < 12; i++)
				{
					k3d::legacy::bicubic_patch* const patch = new k3d::legacy::bicubic_patch();
					patch->material = material;

					for(unsigned long j = 0; j < 16; j++)
						patch->control_points[j] = Mesh.points[TeapotPatches[i][j]-1];

					Mesh.bicubic_patches.push_back(patch);
				}
				// Handle
				for(unsigned long i = 12; i < 16; i++)
				{
					k3d::legacy::bicubic_patch* const patch = new k3d::legacy::bicubic_patch();
					patch->material = material;

					for(unsigned long j = 0; j < 16; j++)
						patch->control_points[j] = Mesh.points[TeapotPatches[i][j]-1];

					Mesh.bicubic_patches.push_back(patch);
				}
				// Spout
				for(unsigned long i = 16; i < 20; i++)
				{
					k3d::legacy::bicubic_patch* const patch = new k3d::legacy::bicubic_patch();
					patch->material = material;

					for(unsigned long j = 0; j < 16; j++)
						patch->control_points[j] = Mesh.points[TeapotPatches[i][j]-1];

					Mesh.bicubic_patches.push_back(patch);
				}
				// Lid
				for(unsigned long i = 20; i < 28; i++)
				{
					k3d::legacy::bicubic_patch* const patch = new k3d::legacy::bicubic_patch();
					patch->material = material;

					for(unsigned long j = 0; j < 16; j++)
						patch->control_points[j] = Mesh.points[TeapotPatches[i][j]-1];

					Mesh.bicubic_patches.push_back(patch);
				}
				// Bottom
				for(unsigned long i = 28; i < 32; i++)
				{
					k3d::legacy::bicubic_patch* const patch = new k3d::legacy::bicubic_patch();
					patch->material = material;

					for(unsigned long j = 0; j < 16; j++)
						patch->control_points[j] = Mesh.points[TeapotPatches[i][j]-1];

					Mesh.bicubic_patches.push_back(patch);
				}

				break;
			}
			case TEACUP :
			{
				// Min : -1.0, max : 1.0
				for(unsigned long i = 0; i < 251; i++)
					Mesh.points.push_back(new k3d::legacy::point(size * k3d::point3(TeacupPoints[i][0], TeacupPoints[i][2], TeacupPoints[i][1])));

				for(unsigned long i = 0; i < 26; i++)
				{
					k3d::legacy::bicubic_patch* const patch = new k3d::legacy::bicubic_patch();
					patch->material = material;

					for(unsigned long j = 0; j < 16; j++)
						patch->control_points[j] = Mesh.points[TeacupPatches[i][j]-1];

					Mesh.bicubic_patches.push_back(patch);
				}

				break;
			}
			case TEASPOON :
			{
				// Min : -1.0, max : 0.221 -> 20/1.221 = 16.38
				for(unsigned long i = 0; i < 256; i++)
					Mesh.points.push_back(new k3d::legacy::point(size * k3d::point3(TeaspoonPoints[i][1], TeaspoonPoints[i][0], TeaspoonPoints[i][2])));

				for(unsigned long i = 0; i < 16; i++)
				{
					k3d::legacy::bicubic_patch* const patch = new k3d::legacy::bicubic_patch();
					patch->material = material;

					for(unsigned long j = 0; j < 16; j++)
						patch->control_points[j] = Mesh.points[TeaspoonPatches[i][j]-1];

					Mesh.bicubic_patches.push_back(patch);
				}

				break;
			}
		}
	}

	void on_update_mesh(k3d::legacy::mesh& Mesh)
	{
	}

	static k3d::iplugin_factory& get_factory()
	{
		static k3d::document_plugin_factory<newell_primitive_implementation, k3d::interface_list<k3d::imesh_source > > factory(
			k3d::uuid(0x274c0cae, 0x2efd5bbf, 0x986a500f, 0xff5e2de6),
			"Newell",
			_("Generates Newell primitives as Bezier patches"),
			"Patches",
			k3d::iplugin_factory::STABLE);

		return factory;
	}

private:
	/// Enumerates supported primitive types
	typedef enum
	{
		TEAPOT,
		TEACUP,
		TEASPOON
	} primitive_t;

	friend std::ostream& operator << (std::ostream& Stream, const primitive_t& Value)
	{
		switch(Value)
		{
			case TEAPOT:
				Stream << "teapot";
				break;
			case TEACUP:
				Stream << "teacup";
				break;
			case TEASPOON:
				Stream << "teaspoon";
				break;
		}
		return Stream;
	}

	friend std::istream& operator >> (std::istream& Stream, primitive_t& Value)
	{
		std::string text;
		Stream >> text;

		if(text == "teapot")
			Value = TEAPOT;
		else if(text == "teacup")
			Value = TEACUP;
		else if(text == "teaspoon")
			Value = TEASPOON;
		else
			k3d::log() << error << k3d_file_reference << ": unknown enumeration [" << text << "]" << std::endl;

		return Stream;
	}

	static const k3d::ienumeration_property::enumeration_values_t& type_values()
	{
		static k3d::ienumeration_property::enumeration_values_t values;
		if(values.empty())
		{
			values.push_back(k3d::ienumeration_property::enumeration_value_t("Teapot", "teapot", "Creates a Newell teapot"));
			values.push_back(k3d::ienumeration_property::enumeration_value_t("Teacup", "teacup", "Creates a Newell teacup"));
			values.push_back(k3d::ienumeration_property::enumeration_value_t("Teaspoon", "teaspoon", "Creates a Newell teaspoon"));
		}

		return values;
	}

	/// Controls the type of primitive to be created
	k3d_data(primitive_t, immutable_name, change_signal, with_undo, local_storage, no_constraint, enumeration_property, with_serialization) m_type;

	/// Generic size control
	k3d_data(double, immutable_name, change_signal, with_undo, local_storage, no_constraint, measurement_property, with_serialization) m_size;
};

/////////////////////////////////////////////////////////////////////////////
// newell_primitive_factory

k3d::iplugin_factory& newell_primitive_factory()
{
	return newell_primitive_implementation::get_factory();
}

} // namespace libk3dprimitives

