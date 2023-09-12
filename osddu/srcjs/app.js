import { OSG } from './OSG'

const App = (props) => {
  return (
    <>
      <OSG
        host={props.host}
        port={props.port}
        width={props.width}
        height={props.height}
        quit={props.quit}
        sampleID={props.sampleID}
      />
    </>
   )
}

export default App
