import sst
import sst.macro

opticalLatency = "0ps"
smallLatency = "100ns"

def makeUniLink(linkType,srcComp,srcId,srcPort,dstComp,dstId,dstPort,outLat=None,inLat=None):
	if not outLat : outLat = inLat
	if not inLat: inLat = outLat
	if not outLat: sys.exit("must specify at least one latency for link")

	linkName = "%s%d:%d->%d:%d" % (linkType,srcId,srcPort,dstId,dstPort)
	link = sst.Link(linkName)
	portName = "output %d %d" % (srcPort, dstPort)
	srcComp.addLink(link,portName,outLat)
	portName = "input %d %d" % (srcPort, dstPort)
	dstComp.addLink(link,portName,inLat)

def makeUniNetworkLink(srcComp,srcId,srcPort,dstComp,dstId,dstPort,outLat=None,inLat=None):
	makeUniLink("network",srcComp,srcId,srcPort,dstComp,dstId,dstPort,outLat,inLat)

def makeBiLink(linkType,comp1,id1,port1,comp2,id2,port2,outLat=None,inLat=None):
	makeUniLink(linkType,comp1,id1,port1,comp2,id2,port2,outLat,inLat)
	makeUniLink(linkType,comp2,id2,port2,comp1,id1,port1,outLat,inLat)

class FlexflyInterconnectSimplified:
	def __init__(self, params):
		self.params = params
		self.simplifiedSwitch = None
		self.system = sst.macro.System(params)
		self.num_nodes = self.system.numNodes()		
		self.nodes = [0]*self.num_nodes #just expands the size of an initial array with 1 element 0 to size 1 * num_nodes
		self.containsOptics = True
		self.optical_topology_params = self.params["topology"]
		self.numGroups = int(self.optical_topology_params["groups"])
		self.switchesPerGroup = int(self.optical_topology_params["switches_per_group"])
		self.numElecSwitches = self.numGroups * self.switchesPerGroup
		self.numSwitches = self.numElecSwitches + 1
		self.nodesPerSwitch = int(self.optical_topology_params["nodes_per_switch"])
		self.optical_switch_radix = self.numGroups
		# really just equals to h in Kim's literature (num of global links per switch)
		self.num_nodes = self.numGroups * self.switchesPerGroup * self.nodesPerSwitch
		self.elecSwitches = []
		self.opticalNetwork = None

	def buildEndpoints(self):
		nodeParams = self.params["node"]
		modelName = nodeParams["model"]
		for i in range(self.num_nodes): 
			self.nodes[i] = sst.Component("Node %d" % i, "macro.%s" % (modelName + "_node"))
			self.nodes[i].addParams(macroToCoreParams(nodeParams))
			self.nodes[i].addParam("id", i)       
		return


	def buildElectricalSwitches(self):
		switchParams = self.params["switch"]
		switchName = "flexfly_electrical" + "_switch"
		#totalSwitchesPerGroup = self.optical_switchesPerGroup + self.switchesPerGroup
		for i in range(self.numElecSwitches):
			#if self.containsOptics and i % (totalSwitchesPerGroup) >= self.switchesPerGroup:
				#continue
			switch = sst.Component("Switch %d" % i, "macro.%s" % switchName)
			switch.addParam("id" , i) 
			switch.addParams(macroToCoreParams(switchParams))
			switch.addParam("switch_type" , "electrical")
			switch.addParam("total_radix", self.switchesPerGroup + self.nodesPerSwitch)
			switch.addParam("switches_per_group", self.switchesPerGroup)
			switch.addParam("num_groups", self.numGroups);
			self.elecSwitches.append(switch) 
		return

	def latency(self, params):
		if params.has_key("latency"):
			return params["latency"]
		elif params.has_key("send_latency"):
			return params["send_latency"]
		else:
			sys.exit("need link latency in parameters")

	def buildTopology(self):
		switchParams = self.params["switch"]
		for i in range(self.numSwitches):
			linkParams = switchParams["link"]
			connections = self.system.switchConnections(i)
			srcSwitch = None
			if i == self.numGroups * self.switchesPerGroup:
				srcSwitch = self.opticalNetwork
			else:
				srcSwitch = self.elecSwitches[i]
			lat = self.latency(linkParams)
			for srcId, dstId, srcOutport, dstInport in connections:
				dstSwitch = None
				assert(srcId <= self.numGroups * self.switchesPerGroup and dstId <= self.numGroups * self.switchesPerGroup)
				if dstId == self.numGroups * self.switchesPerGroup:
					dstSwitch = self.opticalNetwork
				else:
					dstSwitch = self.elecSwitches[dstId]
				makeUniNetworkLink(srcSwitch, srcId, srcOutport, dstSwitch, dstId, dstInport, smallLatency)
		return

	def buildNodeConnections(self):
		for i in range(self.numElecSwitches):
			switch = self.elecSwitches[i]
			for nodeIndex in range(self.nodesPerSwitch):
				index = nodeIndex + (i * self.nodesPerSwitch)
				node = self.nodes[index]
				switchPortIndex = self.switchesPerGroup - 1 + 1 + nodeIndex
				makeUniLink("injection",node,index,0,switch,i,switchPortIndex,smallLatency)
				makeUniLink("ejection",switch,i,switchPortIndex,node,index,0,smallLatency)

	def buildIntraGroupConnections(self):
		inports = [0] * self.numElecSwitches
		outports = [0] * self.numElecSwitches
		for group in range(self.numGroups):
			for i in range(self.switchesPerGroup):
				srcID = i + group * self.switchesPerGroup
				srcSwitch = self.elecSwitches[srcID]
				if i == self.switchesPerGroup - 1:
					break
				for j in range(i + 1, self.switchesPerGroup, 1):
					dstID = j + group * self.switchesPerGroup
					dstSwitch = self.elecSwitches[dstID]
					linkName = "network %d:%d->%d:%d" % (srcID, outports[srcID] , dstID, inports[dstID])
					makeUniLink(linkName, srcSwitch, srcID, outports[srcID], dstSwitch, dstID, inports[dstID], smallLatency)
					linkName = "network %d:%d->%d:%d" % (dstID, outports[dstID] , srcID, inports[srcID])
					makeUniLink(linkName, dstSwitch, dstID, outports[dstID], srcSwitch, srcID, inports[srcID], smallLatency)
					outports[srcID] += 1
					outports[dstID] += 1
					inports[srcID] += 1
					inports[dstID] += 1

	def buildNodeConnections2(self):
		for i in range(self.numGroups * self.switchesPerGroup):
			switch = self.elecSwitches[i]
			for nodeIndex in range(self.nodesPerSwitch):
				index = nodeIndex + (i * self.nodesPerSwitch)
				node = self.nodes[index]
				switchPortIndex = self.switchesPerGroup + nodeIndex
				linkName = "logP %d -> %d" % (index, i)
				link = sst.Link(linkName)
				portName = "in-out %d %d" % (sst.macro.NICLogPInjectionPort, sst.macro.SwitchLogPInjectionPort)
				node.addLink(link, portName, smallLatency)
				portName = "in-out %d %d" % (i, sst.macro.SwitchLogPInjectionPort)
				switch.addLink(link, portName, smallLatency)

	## 
	## Builds the simplified optical network model, in which all the 1 optical link from electrical
	## switches go into 1 single instantiation of the class called optical network
	## This provides an abstraction to prevent us from needing to deal with the intricacies of 
	## routing through the optical network
	##
	def buildOpticalNetwork(self):
		print "Building Optical Network subfunction in Python"
		self.opticalNetwork = sst.Component("flexfly_optical_network %d" % 0, "macro.flexfly_optical_network")
		self.opticalNetwork.addParam("num_groups", self.numGroups)
		self.opticalNetwork.addParam("switches_per_group", self.switchesPerGroup)
		self.opticalNetwork.addParam("num_electrical_switches", self.numElecSwitches)
		self.opticalNetwork.addParam("id", self.numGroups * self.switchesPerGroup)
		return 

	def buildOpticalNetworkLinks(self):
		opticalNetworkPort = 0
		i = 0
		electricalSwitchPort = self.switchesPerGroup - 1
		for elecSwitch in self.elecSwitches:
			linkName = "opticalNetwork %d:%d->%d:%d" % (i,electricalSwitchPort,0,opticalNetworkPort)
			link = sst.Link(linkName)
			makeUniLink("elec->optical", elecSwitch, i, electricalSwitchPort, self.opticalNetwork, 3000, i, smallLatency)
			makeUniLink("elec->optical", self.opticalNetwork, 3000, i, elecSwitch, i, electricalSwitchPort, smallLatency)
			i += 1
			opticalNetworkPort += 1
		return

	def buildLogPNetwork(self):
		import re
		nproc = sst.getMPIRankCount() * sst.getThreadCount()
		switchParams = self.params["switch"]
		linkParams = switchParams["link"]
		ejParams = switchParams["ejection"]
		lat = self.latency(ejParams)
		#gotta multiply the lat by 2
		match = re.compile("(\d+[.]?\d*)(.*)").search(lat)
		if not match:
			sys.exit("improperly formatted latency %s" % lat)
		num, units = match.groups()
		num = eval(num) * 2
		lat = "%8.4f%s" % (num,units.strip())
		switches = []
		for i in range(nproc):
			switch = sst.Component("my_logp %d" % i, "macro.my_logp_switch")
			#switch = sst.Component("logp %d" % i, "macro.logp_switch")
			switch.addParams(macroToCoreParams(switchParams))
			switch.addParam("id", i)
			switches.append(switch)
			#NEWLY ADDED BEGIN
			switch.addParam("electrical_bandwidth" , "10Gb/s")
			switch.addParam("optical_bandwidth", "1Gb/s")
			#NEWLY ADDED END
		for i in range(nproc):
			sw_i = switches[i]
			for j in range(nproc):
				sw_j = switches[j]
				if i==j: continue
				linkName = "logPnetwork%d->%d" % (i,j)
				link = sst.Link(linkName)
				portName = "in-out %d %d" % (j, sst.macro.SwitchLogPNetworkPort)
				sw_i.addLink(link, portName, lat)
				portName = "in-out %d %d" % (i, sst.macro.SwitchLogPNetworkPort)
				sw_j.addLink(link, portName, lat)			
		
		for i in range(self.num_nodes):
			injSW = self.system.nodeToLogPSwitch(i)
			ep = self.nodes[i]
			sw = switches[injSW]
			linkName = "logPinjection%d->%d" % (i, injSW)
			link = sst.Link(linkName)
			portName = "in-out %d %d" % (sst.macro.NICLogPInjectionPort, sst.macro.SwitchLogPInjectionPort)
			#portName = "injection %d %d" % (sst.macro.NICLogPInjectionPort, sst.macro.SwitchLogPInjectionPort)
			ep.addLink(link, portName, smallLatency) #put no latency here
			portName = "in-out %d %d" % (i, sst.macro.SwitchLogPInjectionPort)
			#portName = "ejection %d %d" % (i, sst.macro.SwitchLogPInjectionPort)
			sw.addLink(link, portName, smallLatency)
	

	def build(self, islogP):
		self.buildEndpoints()
		if not islogP:
			self.buildElectricalSwitches()
			#self.buildIntraGroupConnections()
			self.buildOpticalNetwork()
			#self.buildOpticalNetworkLinks()
			self.buildTopology()
			self.buildNodeConnections()
			self.buildLogPNetwork()			
		else:
			self.buildEndpoints()
			self.buildSimplifiedTopology(smallLatency)
			self.buildLogPNetwork()
		
		return
		

## Returns the command line argv in terms of a vector that is 0-indexed
def readCmdLineParams():
	import sys
	return sst.macro.readParams(sys.argv)

def redoSubParams_impl(nsArr, theDict, allParams):
	for key in theDict:
		val = theDict[key]
		if isinstance(val, dict):
			newNsArr = nsArr[:]
			newNsArr.append(key)
			redoSubParams_impl(newNsArr, val, allParams)
		else:
			paramArr = nsArr[:]
			paramArr.append(key)
			newParam = ".".join(paramArr)
			allParams.append((newParam, val))

def macroToCoreParams(theDict):
	allParams = []
	redoSubParams_impl([], theDict, allParams)
	newDict = {}
	for key, val in allParams:
		newDict[key] = val
	return newDict

def setupTopologySimplified():
	import sys
	params = readCmdLineParams()
	nodeParams = params["node"]
	swParams = params["switch"]
	
	builtinApps = [
   		"apitest",
       	"global_test",
       	"hello_world",
       	"mpi_coverage",
       	"mpi_ping_all",
       	"mpi_print_nodes",
       	"mpi_topology",
       	"parsedumpi",
       	"sstmac_mpi_testall",
       	"traffic_matrix",
       	"user_app_cxx_empty_main",
       	"user_app_cxx_full_main",
	]

	for appIdx in range(10):
		appKey = "app%d" % (appIdx)
		#nodeParams = params["node"]
		if params.has_key(appKey):
			appParams = params[appKey]
			nodeParams[appKey] = appParams
			appName = appParams["name"]
			if not appName in builtinApps:
				cmd = "import sst.%s" % appName
				exec(cmd)
			del params[appKey]

	debugList = []
	if params.has_key("debug"):
		debugList = params["debug"].strip().split()
	for i in range(len(sys.argv)):
		if sys.argv[i] == "-d" or sys.argv[i] == "--debug":
			debugList.extend(sys.argv[i+1].split(","))

	icParams = {}
	icParams["topology"] = params["topology"]
	nodeParams["interconnect"] = icParams
	if debugList:
		nodeParams["debug"] = " ".join(debugList)
	swParams["topology"] = params["topology"]

	swParams["topology"] = "torus"

	#move every param in the global namespace 
	#into the individal namespaces
	for ns in "node", "switch":
		nsParams = params[ns]
		for key in params:
			val = params[key]
			if isinstance(val, str):
				if not nsParams.has_key(key):
					nsParams[key] = val    
	interconnect = FlexflyInterconnectSimplified(params)
	isLogP = False
	interconnect.build(isLogP)
