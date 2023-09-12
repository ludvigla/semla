import React, { useEffect, useState } from "react";
import OpenSeaDragon from "openseadragon";
import "./opensea-dragon-svg-overlay";
import * as d3 from "d3";
import useKeyPress from "./UseKeyPress";

window.shouldQuit = false

const OSG = (props) => {

  // Create state variables for viewer, overlay and d3 nodes (svg)
  const [imageInfo, setImageInfo] = useState(null)
  const [viewer, setViewer] = useState(null)
  const [overlay, setOverlay] = useState(null)
  const [d3Node, setd3Node] = useState(null)
  const [dataset, setDataset] = useState(null)

  const getImageInfo = () => {
    fetch(`http://${props.host}:${props.port}/image_info_${props.sampleID}.json`)
      .then(function (response) {
        return response.json();
      })
      .then(function (myJson) {
        setImageInfo(myJson);
      });
  };

  // load image info and data on render
  useEffect(() => {
      getImageInfo()
      getData()
  }, [])

  const getData = () => {
    fetch(`http://${props.host}:${props.port}/network_Visium_${props.sampleID}.json`)
      .then(function (response) {
        return response.json();
      })
      .then(function (myJson) {
        setDataset(myJson);
      });
  };

 // Function to initiate the OpenSeaDragon viewer
 const InitOpenseadragon = () => {
  viewer && viewer.destroy();
  setViewer(
    OpenSeaDragon({
      id: "openseadragon",
      tileSources: {
        height: imageInfo.image_height,
        width: imageInfo.image_width,
        tileSize: imageInfo.tilesize,
        minLevel: imageInfo.minZoomLevel,
        maxLevel: imageInfo.maxZoomLevel,
        getTileUrl: (level, x, y) => {
          return `http://${props.host}:${props.port}/tiles${props.sampleID}/${level}/${x}_${y}.jpg`;
        }
      },
      defaultZoomLevel: imageInfo.minZoomLevel,
      minZoomLevel: imageInfo.minZoomLevel,
      maxZoomLevel: imageInfo.maxZoomLevel + 2,
      visibilityRatio: 1.0,
      minZoomImageRatio: 1.0,
      showNavigationControl: false
    })
  );
 };

 // Initiate a new osd viewer once image info has been loaded
 useEffect(() => {
  if (imageInfo) {
    InitOpenseadragon();
  }
  return () => {
    viewer && viewer.destroy();
  };
  // eslint-disable-next-line react-hooks/exhaustive-deps
}, [imageInfo]);

// Add a svg overlay to the osd viewer after the viewer is initialized
useEffect(() => {
  if (viewer) {
    setOverlay(viewer.svgOverlay());
  }
}, [viewer]);


// Create d3 nodes for the overlay
useEffect(() => {
  if (overlay) {
    setd3Node(d3.select(overlay.node()));
  }
}, [overlay]);

// Draw spots on the svg overlay
useEffect(() => {
  if (dataset) {
    //console.log(dataset);

    // Draw links between spots
    d3Node
      .selectAll("line")
      .data(dataset.links)
      .enter()
      .append('svg:line')
      .attr("x1", function (d) {
        return d.x;
      })
      .attr("y1", function (d) {
        return d.y;
      })
      .attr("x2", function (d) {
        return d.x_end;
      })
      .attr("y2", function (d) {
        return d.y_end;
      })
      .attr("shape-rendering", "optimizeSpeed")
      .attr("style", function(d) {
        if (d.keep) {
          return "stroke: black"
        } else {
          return "stroke: red"
        }
      })
      .style('stroke-width', 0.001)
      .on("mouseover", function (e) {
        if (e.shiftKey) {
          const obj = e.srcElement.__data__;
          let copy = dataset;
          copy.links[obj.index - 1].keep = false;
          setDataset(copy);
          d3.select(this)
            .style("stroke", "red")
        }
        if (e.ctrlKey) {
          const obj = e.srcElement.__data__;
          let copy = dataset;
          copy.links[obj.index - 1].keep = true;
          setDataset(copy);
          d3.select(this)
            .style("stroke", "black")
            .attr("remove", false)
        }
      })

    // draw circles for each spot
    d3Node
      .selectAll("circle")
      .data(dataset.nodes)
      .enter()
      .append("g")
      .append("circle")
      .attr("cx", function (d) {
        return d.x;
      })
      .attr("cy", function (d) {
        return d.y;
      })
      .attr("shape-rendering", "optimizeSpeed")
      .attr("r", 0.003)
      .style("fill", "white")
      .attr("opacity", function (d) {
        return 0.8;
      })

  }
}, [dataset, d3Node]);

 // create a key press
 const shiftPress = useKeyPress("Shift");
 const ctrlPress = useKeyPress("Control");

 useEffect(() => {
  if (shiftPress) {
    document.getElementById("openseadragon").style.cursor = "url('data:image/x-icon;base64,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'),auto"
  } else if (ctrlPress) {
    document.getElementById("openseadragon").style.cursor = "url('data:image/x-icon;base64,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'), auto";
  } else {
    document.getElementById("openseadragon").style.cursor = "default";
  }
 }, [shiftPress, ctrlPress]);

  // Listen to changes in window.shouldQuit
  useEffect(() => {
    if (props.quit) {
      window.savedData = dataset
    }
  }, [props.quit, dataset])

  return (
    <>
      <div
        id="openseadragon"
        style={{
          height: props.height,
          width: props.width,
        }}
      ></div>
    </>
  )
};

export { OSG };
