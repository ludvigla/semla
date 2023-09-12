import React, { useEffect, useState } from "react"
import OpenSeaDragon from "openseadragon"
import "./opensea-dragon-svg-overlay"
import * as d3 from "d3"
import "./OSG.css"
import Lasso from "./components/Lasso"
import Slider from "./components/Slider"
import utils from "./utils/utils"
import d3_utils from "./utils/d3_utils"

const OSG = (props) => {
  // Create state variable for the current selected path (lasso tool)
  // and pass data from the canvas to the OSD viewer
  const [curPath, setCurPath] = useState(false)
  const passData = (data) => {
    setCurPath(data)
  }

  // Create state variable for the current selected path to erase (lasso tool)
  // and pass data from the canvas to the OSD viewer
  const [curErasePath, setCurErasePath] = useState(false)
  const passEraseData = (data) => {
    setCurErasePath(data)
  }

  // Create state variables for viewer, overlay and d3 nodes (svg)
  const [imageInfo, setImageInfo] = useState(null)
  const [viewer, setViewer] = useState(null)
  const [overlay, setOverlay] = useState(null)
  const [d3Node, setd3Node] = useState(null)
  const [dataset, setDataset] = useState(null)
  // const [opacity, setOpacity] = useState(1)

  // load image info and data on render
  useEffect(() => {
    utils.getImageInfo(props.host, props.port, setImageInfo, props.sampleID)
    utils.getData(props.host, props.port, props.values, props.opacities, props.sampleID, setDataset)
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [props.sampleID])

  // Initiate a new osd viewer once image info has been loaded
  useEffect(() => {
    if (imageInfo) {
      utils.InitOpenseadragon(
        props.host,
        props.port,
        props.sampleID,
        viewer,
        setViewer,
        OpenSeaDragon,
        imageInfo
      )
    }
    return () => {
      viewer && viewer.destroy()
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [imageInfo])

  // Add a svg overlay to the osd viewer after the viewer is initialized
  useEffect(() => {
    if (viewer) {
      setOverlay(viewer.svgOverlay())
    }
  }, [viewer, props.sampleID])

  // Create d3 nodes for the overlay
  useEffect(() => {
    if (overlay && dataset) {
      setd3Node(d3.select(overlay.node()))
    }
  }, [overlay, dataset, props.sampleID])

  // Draw spots on the svg overlay
  useEffect(() => {
    if (d3Node && dataset) {
      d3_utils.draw_circles(d3Node, dataset)
    }
  }, [d3Node, dataset, props.sampleID])

  // Change opacity
  useEffect(() => {
    if (d3Node && dataset) {
      d3_utils.change_opacity(d3Node, dataset, props.opacity, props.scaleByOpacity, props.opacities)
    }
  }, [d3Node, props.opacity, dataset, props.scaleByOpacity, props.opacities])

  // Change colors
  useEffect(() => {
    if (d3Node && dataset) {
      if (props.isNumeric) {
        d3_utils.color_by_numeric(d3Node, dataset, props.colors, props.values, props.range)
      } else {
       d3_utils.color_by_categorical(d3Node, dataset, props.levels, props.colors, props.values)
      }
    }
  }, [d3Node, dataset, props.isNumeric, props.levels, props.colors, props.values, props.range])

  // Select spots
  useEffect(() => {
    if (curPath) {
      var prps = utils.getViewerProps(viewer)
      let datacopy = dataset
      datacopy = d3_utils.isPointInsidePolygon(datacopy, prps, curPath, true, props.width, props.height)
      setDataset(datacopy)
      if (d3Node && dataset) {
        d3_utils.change_stroke_selected(d3Node, dataset)
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [curPath, d3Node, dataset, viewer])

  // Deselect spots
  useEffect(() => {
    if (curErasePath && dataset) {
      var prps = utils.getViewerProps(viewer)
      let datacopy = dataset
      datacopy = d3_utils.isPointInsidePolygon(datacopy, prps, curErasePath, false, props.width, props.height)
      setDataset(datacopy)
      if (d3Node && dataset) {
        d3_utils.change_stroke_deselected(d3Node, dataset)
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [curErasePath, dataset, viewer])

  // Clear lasso select
  useEffect(() => {
    if (d3Node && dataset) {
      if (!props.useLasso) {
        setCurPath(false)
        setCurErasePath(false)
        d3Node
          .selectAll('circle')
          .data(dataset.nodes)
          .style('stroke-width', (d) => 0)
        const dataselection = d3Node
          .selectAll('circle')
          .data(dataset.nodes)
          .filter((d) => d.selected)
          .data()
        const barcodes = []
        dataselection.forEach((user) => barcodes.push(user.barcode))
        window.curSelection = barcodes
        const datacopy = dataset
        datacopy.nodes = datacopy.nodes.map((node) => {
          node.selected = false
          return node
        })
        setDataset(datacopy)
      }
    }
  }, [d3Node, dataset, props.useLasso])

  return (
    <>
      <div
        className="container"
        style={{
          height: `${props.width}px`,
          width: `${props.height}px`
        }}
      >
        <div
          id="openseadragon"
          className="content"
          style={{
            height: `${props.width}px`,
            width: `${props.height}px`,
          }}
        ></div>
        <div className="overlay">
          {props.useLasso ? (
              <div>
                <Lasso
                  passData={passData}
                  passEraseData={passEraseData}
                  width={`${props.width}px`}
                  height={`${props.height}px`}
                />
              </div>
            ) : null}
        </div>
      </div>
    </>
  )
};

export { OSG };
