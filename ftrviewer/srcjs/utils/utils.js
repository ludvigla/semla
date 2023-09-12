const getImageInfo = (host, port, setImageInfo, sampleID) => {
  fetch(`http://${host}:${port}/image_info_${sampleID}.json`)
    .then(function (response) {
      return response.json()
    })
    .then(function (myJson) {
      setImageInfo(myJson)
    })
}

const getData = (host, port, values, opacities, sampleID, setDataset) => {
  fetch(
    `http://${host}:${port}/coords_Visium_${sampleID}.json`
  )
    .then(function (response) {
      return response.json()
    })
    .then(function (myJson) {
      myJson.nodes = myJson.nodes.map((node, index) => {
        node.value = values[index]
        node.opacity = opacities[index]
        return node
      })
      setDataset(myJson)
    })
}

// Function to initiate the OpenSeaDragon viewer
const InitOpenseadragon = (
  host,
  port,
  sampleID,
  viewer,
  setViewer,
  OpenSeaDragon,
  imageInfo
) => {
  viewer && viewer.destroy()
  setViewer(
    OpenSeaDragon({
      id: 'openseadragon',
      tileSources: {
        height: imageInfo.image_height,
        width: imageInfo.image_width,
        tileSize: imageInfo.tilesize,
        minLevel: imageInfo.minZoomLevel,
        maxLevel: imageInfo.maxZoomLevel,
        getTileUrl: (level, x, y) => {
          return `http://${host}:${port}/tiles${sampleID}/${level}/${x}_${y}.jpg`
        },
      },
      defaultZoomLevel: imageInfo.minZoomLevel,
      minZoomLevel: imageInfo.minZoomLevel,
      maxZoomLevel: imageInfo.maxZoomLevel + 2,
      visibilityRatio: 1.0,
      minZoomImageRatio: 1.0,
      showNavigationControl: false,
    })
  )
}

const getViewerProps = (viewer) => {
  var curBounds = viewer.viewport.getBounds()
  var prps = {}
  prps.x = curBounds.x
  prps.y = curBounds.y
  prps.width = curBounds.width
  prps.height = curBounds.height
  return prps
}

const convertRange = (value, r1, r2) => {
  return ((value - r1[0]) * (r2[1] - r2[0])) / (r1[1] - r1[0]) + r2[0]
}

export default { getImageInfo, getData, InitOpenseadragon, getViewerProps, convertRange }
