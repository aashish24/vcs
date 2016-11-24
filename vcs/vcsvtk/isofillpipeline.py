from .pipeline2d import Pipeline2D
import fillareautils

import numpy
import vcs
import vtk


def vtklutTomplColor(vtklut):
    from matplotlib.colors import LinearSegmentedColormap

    vtklut.ForceBuild()
    vtkcolors = vtklut.GetTable()
    colors = []
    print 'SetNumberOfTableValues ', vtklut.GetNumberOfTableValues()
    noOfColors = vtklut.GetNumberOfTableValues()

    for i in range(0, noOfColors):
        tupl = vtklut.GetTableValue(i)
        print tupl
        colors.append(tupl)
    cmap = LinearSegmentedColormap.from_list('mycmap', colors, N=noOfColors)
    return cmap

def writePdf(actors):
    import vtk
    from numpy import zeros
    import matplotlib.pyplot as plt

    for actorList in actors:
        actor = actorList[0]
        data = actor.GetMapper().GetInput()

        filter = vtk.vtkTriangleFilter()
        filter.SetInputData(data)
        filter.Update()

        data = filter.GetOutput()

        triangles = data.GetPolys()
        points = data.GetPoints()

        ntri = triangles.GetNumberOfCells()
        npts = points.GetNumberOfPoints()

        tri = zeros((ntri, 3))
        x = zeros(npts)
        y = zeros(npts)
        z = zeros(npts)

        # ux = zeros(nvls)
        # uy = zeros(nvls)

        idlist = vtk.vtkIdList()

        for i in xrange(0, ntri):
            cell = triangles.GetNextCell(idlist)
            tri[i, 0] = idlist.GetId(0)
            tri[i, 1] = idlist.GetId(1)
            tri[i, 2] = idlist.GetId(2)

        for i in xrange(npts):
            pt = points.GetPoint(i)
            x[i] = pt[0]
            y[i] = pt[1]
            z[i] = pt[2]

        # for i in xrange(0, nvls):
        #     U = vels.GetTuple(i)
        #     ux[i] = U[0]
        #     uy[i] = U[1]

        # Mesh
        # plt.figure(figsize=(8, 8))
        # plt.triplot(x, y, tri)
        # plt.gca().set_aspect('equal')
        # plt.show()


        mapper = vtk.vtkCellDataToPointData()
        mapper.AddInputData(data)
        mapper.Update()
        data = mapper.GetOutput()
        vels = data.GetPointData().GetArray(2)
        nvls = vels.GetNumberOfTuples()

        ux = zeros(nvls)
        uy = zeros(nvls)

        print vels.GetNumberOfTuples()

        for i in xrange(0, nvls):
            U = vels.GetTuple(i)
            ux[i] = U[0]
            # print ux[i]


        from matplotlib.backends.backend_pdf import PdfPages

        # with PdfPages('multipage_pdf.pdf') as pdf:
        plt.figure()
        cmap = vtklutTomplColor(actor.GetMapper().GetLookupTable())
        plt.tricontourf(x, y, tri, ux, 4, cmap=cmap, antialiased=True)
        plt.show()
            # pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(data)
    writer.SetFileName("foo.vtk")
    writer.Write()


class IsofillPipeline(Pipeline2D):

    """Implementation of the Pipeline interface for VCS isofill plots."""

    def __init__(self, gm, context_):
        super(IsofillPipeline, self).__init__(gm, context_)
        self._needsCellData = False

    def _updateContourLevelsAndColors(self):
        self._updateContourLevelsAndColorsGeneric()

    def _plotInternal(self):
        """Overrides baseclass implementation."""

        preppedCountours = self._prepContours()
        tmpLevels = preppedCountours["tmpLevels"]
        tmpIndices = preppedCountours["tmpIndices"]
        tmpColors = preppedCountours["tmpColors"]
        tmpOpacities = preppedCountours["tmpOpacities"]
        style = self._gm.fillareastyle

        luts = []
        cots = []
        mappers = []
        _colorMap = self.getColorMap()

        plotting_dataset_bounds = self.getPlottingBounds()
        x1, x2, y1, y2 = plotting_dataset_bounds

        for i, l in enumerate(tmpLevels):
            # Ok here we are trying to group together levels can be, a join
            # will happen if: next set of levels continues where one left off
            # AND pattern is identical
            mapper = vtk.vtkPolyDataMapper()
            lut = vtk.vtkLookupTable()
            cot = vtk.vtkBandedPolyDataContourFilter()
            cot.ClippingOn()
            cot.SetInputData(self._vtkPolyDataFilter.GetOutput())
            cot.SetNumberOfContours(len(l))
            cot.SetClipTolerance(0.)
            for j, v in enumerate(l):
                cot.SetValue(j, v)
            cot.Update()

            cots.append(cot)
            mapper.SetInputConnection(cot.GetOutputPort())
            lut.SetNumberOfTableValues(len(tmpColors[i]))
            for j, color in enumerate(tmpColors[i]):
                r, g, b, a = self.getColorIndexOrRGBA(_colorMap, color)
                if style == 'solid':
                    tmpOpacity = tmpOpacities[j]
                    if tmpOpacity is None:
                        tmpOpacity = a / 100.
                    else:
                        tmpOpacity = tmpOpacities[j] / 100.
                    lut.SetTableValue(j, r / 100., g / 100., b / 100., tmpOpacity)
                else:
                    lut.SetTableValue(j, 1., 1., 1., 0.)
            luts.append([lut, [0, len(l) - 1, True]])
            mapper.SetLookupTable(lut)
            minRange = 0
            maxRange = len(l) - 1
            if (i == 0 and self._scalarRange[0] < l[0]):
                # band 0 is from self._scalarRange[0] to l[0]
                # we don't show band 0
                minRange += 1
            mapper.SetScalarRange(minRange, maxRange)
            mapper.SetScalarModeToUseCellData()
            mappers.append(mapper)

        self._resultDict["vtk_backend_luts"] = luts
        if len(cots) > 0:
            self._resultDict["vtk_backend_contours"] = cots

        numLevels = len(self._contourLevels)
        if mappers == []:  # ok didn't need to have special banded contours
            mapper = vtk.vtkPolyDataMapper()
            mappers = [mapper]
            # Colortable bit
            # make sure length match
            while len(self._contourColors) < len(self._contourLevels):
                self._contourColors.append(self._contourColors[-1])

            lut = vtk.vtkLookupTable()
            lut.SetNumberOfTableValues(numLevels)
            for i in range(numLevels):
                r, g, b, a = self.getColorIndexOrRGBA(_colorMap, self._contourColors[i])
                lut.SetTableValue(i, r / 100., g / 100., b / 100., a / 100.)
            mapper.SetLookupTable(lut)
            if numpy.allclose(self._contourLevels[0], -1.e20):
                lmn = self._min - 1.
            else:
                lmn = self._contourLevels[0]
            if numpy.allclose(self._contourLevels[-1], 1.e20):
                lmx = self._max + 1.
            else:
                lmx = self._contourLevels[-1]
            mapper.SetScalarRange(lmn, lmx)
            self._resultDict["vtk_backend_luts"] = [[lut, [lmn, lmx, True]]]

        if self._maskedDataMapper is not None:
            mappers.insert(0, self._maskedDataMapper)

        # And now we need actors to actually render this thing
        actors = []
        patternActors = []
        ct = 0
        vp = self._resultDict.get('ratio_autot_viewport',
                                  [self._template.data.x1, self._template.data.x2,
                                   self._template.data.y1, self._template.data.y2])
        dataset_renderer = None
        xScale, yScale = (1, 1)
        for mapper in mappers:
            act = vtk.vtkActor()
            act.SetMapper(mapper)

            patact = None
            # TODO see comment in boxfill.
            if mapper is self._maskedDataMapper:
                actors.append([act, self._maskedDataMapper, plotting_dataset_bounds])
            else:
                actors.append([act, plotting_dataset_bounds])

                # Since pattern creation requires a single color, assuming the first
                c = self.getColorIndexOrRGBA(_colorMap, tmpColors[ct][0])

                # The isofill actor is scaled by the camera, so we need to use this size
                # instead of window size for scaling the pattern.
                viewsize = (x2 - x1, y2 - y1)
                patact = fillareautils.make_patterned_polydata(mapper.GetInput(),
                                                               fillareastyle=style,
                                                               fillareaindex=tmpIndices[ct],
                                                               fillareacolors=c,
                                                               fillareaopacity=tmpOpacities[ct],
                                                               size=viewsize)

                if patact is not None:
                    patternActors.append(patact)

                # increment the count
                ct += 1

            # create a new renderer for this mapper
            # (we need one for each mapper because of cmaera flips)
            dataset_renderer, xScale, yScale = self._context().fitToViewport(
                act, vp,
                wc=plotting_dataset_bounds, geoBounds=self._vtkDataSet.GetBounds(),
                geo=self._vtkGeoTransform,
                priority=self._template.data.priority,
                create_renderer=(mapper is self._maskedDataMapper or dataset_renderer is None))
        for act in patternActors:
            self._context().fitToViewport(
                act, vp,
                wc=plotting_dataset_bounds, geoBounds=self._vtkDataSet.GetBounds(),
                geo=self._vtkGeoTransform,
                priority=self._template.data.priority,
                create_renderer=True)
            actors.append([act, plotting_dataset_bounds])

        self._resultDict["vtk_backend_actors"] = actors

        writePdf(actors)

        t = self._originalData1.getTime()
        if self._originalData1.ndim > 2:
            z = self._originalData1.getAxis(-3)
        else:
            z = None
        kwargs = {"vtk_backend_grid": self._vtkDataSet,
                  "dataset_bounds": self._vtkDataSetBounds,
                  "plotting_dataset_bounds": plotting_dataset_bounds,
                  "vtk_backend_geo": self._vtkGeoTransform}
        if ("ratio_autot_viewport" in self._resultDict):
            kwargs["ratio_autot_viewport"] = vp
        self._resultDict.update(self._context().renderTemplate(
            self._template,
            self._data1,
            self._gm, t, z, **kwargs))
        legend = getattr(self._gm, "legend", None)

        if self._gm.ext_1:
            if isinstance(self._contourLevels[0], list):
                if numpy.less(abs(self._contourLevels[0][0]), 1.e20):
                    # Ok we need to add the ext levels
                    self._contourLevels.insert(
                        0, [-1.e20, self._contourLevels[0][0]])
            else:
                if numpy.less(abs(self._contourLevels[0]), 1.e20):
                    # need to add an ext
                    self._contourLevels.insert(0, -1.e20)
        if self._gm.ext_2:
            if isinstance(self._contourLevels[-1], list):
                if numpy.less(abs(self._contourLevels[-1][1]), 1.e20):
                    # need ext
                    self._contourLevels.append([self._contourLevels[-1][1],
                                                1.e20])
            else:
                if numpy.less(abs(self._contourLevels[-1]), 1.e20):
                    # need exts
                    self._contourLevels.append(1.e20)

        self._resultDict.update(
            self._context().renderColorBar(self._template, self._contourLevels,
                                           self._contourColors, legend,
                                           self.getColorMap(),
                                           style=style,
                                           index=self._gm.fillareaindices,
                                           opacity=self._gm.fillareaopacity))

        if self._context().canvas._continents is None:
            self._useContinents = False
        if self._useContinents:
            projection = vcs.elements["projection"][self._gm.projection]
            continents_renderer, xScale, yScale = self._context().plotContinents(
                plotting_dataset_bounds, projection,
                self._dataWrapModulo,
                vp, self._template.data.priority,
                vtk_backend_grid=self._vtkDataSet,
                dataset_bounds=self._vtkDataSetBounds)
