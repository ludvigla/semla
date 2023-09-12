import * as d3 from "d3"
import utils from "./utils"
import Paper from "paper"

// draw circles for each spot
const draw_circles = (d3Node, dataset) => {
  d3Node
    .selectAll('circle')
    .data(dataset.nodes)
    .enter()
    .append('g')
    .append('circle')
    .attr('cx', function (d) {
      return d.x
    })
    .attr('cy', function (d) {
      return d.y
    })
    .attr('r', 0.003)
}

// change opacity of spots
const change_opacity = (d3Node, dataset, opacity, scaleByOpacity, opacities) => {
  !scaleByOpacity
    ? d3Node
        .selectAll('circle')
        .data(dataset.nodes)
        .attr('opacity', function (d) {
          return opacity
        })
    : d3Node
        .selectAll('circle')
        .data(dataset.nodes)
        .attr('opacity', function (d, i) {
          return opacities[i]
        })
}

// Color spots by numeric values
const color_by_numeric = (d3Node, dataset, colors, values, valuesrange) => {
  const color_values = d3.scaleQuantize().domain(valuesrange).range(colors)

  d3Node
    .selectAll('circle')
    .data(dataset.nodes)
    .style('fill', (d, i) => color_values(values[i]))
}

// Color spots by categorical values
const color_by_categorical = (d3Node, dataset, levels, colors, values) => {
  const colors_values = d3
          .scaleOrdinal()
          .domain(levels)
          .range(colors)

  d3Node
    .selectAll('circle')
    .data(dataset.nodes)
    .style('fill', (d, i) => colors_values(values[i]))
}

// Change stroke for lasso
const change_stroke_selected = (d3Node, dataset) => {
  d3Node
    .selectAll('circle')
    .data(dataset.nodes)
    .style('stroke-width', (d) => (d.selected ? 0.002 : 0))
    .style('stroke', '#4682B4')
}

// Change stroke deselected
const change_stroke_deselected = (d3Node, dataset) => {
  d3Node
    .selectAll('circle')
    .data(dataset.nodes)
    .style('stroke-width', (d) => (d.selected ? 0.002 : 0))
}

// Check if a point is inside a polygon
const isPointInsidePolygon = (datacopy, prps, selectedPath, isSelected, width, height) => {
  datacopy.nodes = datacopy.nodes.map((d) => {
    var newPoint = new Paper.Point(
      utils.convertRange(
        d.x,
        [prps.x, prps.width + prps.x],
        [0, width]
      ),
      utils.convertRange(
        d.y,
        [prps.y, prps.height + prps.y],
        [0, height]
      )
    )
    d.selected = selectedPath.contains(newPoint) ? isSelected : d.selected
    return d
  })
  return datacopy
}

const d3_utils = {
  draw_circles,
  change_opacity,
  color_by_numeric,
  color_by_categorical,
  change_stroke_selected,
  change_stroke_deselected,
  isPointInsidePolygon
}

export default d3_utils
