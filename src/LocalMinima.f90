subroutine LocalMinima (stack,rec,bc,ndon,donor,h,length,nx,ny,dx,dy)

  ! subroutine to compute and remove local inima by recomputing the receiver connectivty
  ! using Guillaume Cordonnier's algorithm as helped by Benoit Bovy and debuged with Jean Braun
  ! This is (at first) a fortran translation of a python routine provided to me by Gullaume
  ! on October 24 2017

  ! in input stack(nx*ny) : stack order as defined by Braun and Willett (2013)
  !          rec(nx*ny) : receiver node information as defined by Braun and Willet (2013)
  !          ndon(nx*ny) : number of donors per node
  !          donor(8,nx*ny) : list of donors per node
  !          h(nx*ny) : nodal heigh field
  !          length(nx*ny) : distance between node and receiver (per node)
  !          nx-ny : dimension of rectangular mesh in x- and y-directions
  !          dx-dy : grid spacing in the x- and y-directions

  ! in output stack : updated stack
  !           rec : updated receiver information
  !           ndon : updated number of donor information
  !           donor : updated donor list
  !           length : updated distance between donot and receiver

  ! The updated receiver list takes into account the existence of local minima in the original receiver list,
  ! i.e. nodes taht are their own receivers. The algorithm finds the geometry of lakes centered on these local
  ! minima as well as their sills. By connecting the lakes through their sills, the algorithm finds the steepest
  ! path (for all nodes) to a node on the boundary. A fake height arrays is created such that lakes
  ! are drained in the order given by the new receiver information. The algorithm is O(n) + (NlogN) ou N is the
  ! number of lakes

  ! Note that the user can use the routines loc_min_3__find_receivers, loc_min_3_find_donors and loc_min_3_find_stack
  ! provided below to compute the receiver, donor information and the stack before calling this routine.

  implicit none

  integer nx,ny
  integer stack(nx*ny),rec(nx*ny),ndon(nx*ny),donor(8,nx*ny)
  double precision h(nx*ny),length(nx*ny),dx,dy
  logical bc(nx*ny)

  integer, dimension(:), allocatable :: basins,outlets,tree,basin_stack
  integer, dimension(:,:), allocatable :: conn_basins,conn_nodes,sills,junki
  double precision, dimension(:), allocatable :: conn_weights,junkd
  logical, dimension(:), allocatable :: active_nodes
  logical continuous_flow,continuous_flow_v2
  integer nbasins,basin0,nconn,nconn_max,tree_size

  integer nn,k,nlocmin

  !continuous_flow = .true.
  continuous_flow = .false.
  continuous_flow_v2 = .true.

  nn=nx*ny

  allocate(basins(nn),outlets(nn))

  call compute_basins (stack,rec,basins,outlets,nn,nbasins)

  nconn_max=nbasins*6
  allocate (conn_basins(nconn_max,2),conn_nodes(nconn_max,2),conn_weights(nconn_max))
  allocate(active_nodes(nn))
  active_nodes=.not.bc

  nlocmin=0
  do k=1,nn
    if (rec(k).eq.k .and. active_nodes(k)) nlocmin=nlocmin+1
  enddo
  !if (nlocmin.gt.0) print*,'nlocmin before',nlocmin

  call connect_basins (nbasins,basins,outlets,rec,stack,active_nodes,h,nx,ny, &
  basin0,nconn,conn_basins,conn_nodes,conn_weights,nconn_max)

  allocate (junki(nconn_max,2))
  junki=conn_basins
  deallocate(conn_basins)
  allocate(conn_basins(nconn,2))
  conn_basins(1:nconn,:)=junki(1:nconn,:)
  junki=conn_nodes
  deallocate(conn_nodes)
  allocate(conn_nodes(nconn,2))
  conn_nodes(1:nconn,:)=junki(1:nconn,:)
  deallocate (junki)

  allocate (junkd(nconn_max))
  junkd=conn_weights
  deallocate(conn_weights)
  allocate(conn_weights(nconn))
  conn_weights(1:nconn)=junkd(1:nconn)
  deallocate (junkd)

  allocate (tree(nbasins-1))
  tree_size=0

  call mst_kruskal(conn_weights,conn_basins,nbasins,nconn,tree,tree_size)

  allocate (sills(nbasins,2))
  allocate (basin_stack(nbasins))

  call order_tree (conn_basins(1:nconn,1),conn_basins(1:nconn,2),conn_nodes(1:nconn,1),conn_nodes(1:nconn,2), &
  nconn,nbasins,basin0,tree,tree_size,sills,basin_stack,continuous_flow)

  if (.not.continuous_flow) then
    if (continuous_flow_v2) then
      call correct_receivers_v2 (rec,length,outlets,conn_basins,conn_nodes,tree,h,nx,ny,dx,dy, &
      nbasins,nconn,tree_size)
    else
      call correct_receivers (rec,length,outlets,conn_basins,conn_nodes,tree,h,nx,ny,dx,dy, &
      nbasins,nconn,tree_size)
    endif
  else
    allocate (junkd(nn))
    junkd = h
    call update_fake_topography (sills,basin_stack,basins,h,nx,ny,dx,dy,nbasins)
    call loc_min_3_find_receivers (h,rec,length,bc,nx,ny,dx,dy)
    h = junkd
    deallocate (junkd)
  endif

  nlocmin=0
  do k=1,nn
    if (rec(k).eq.k.and.active_nodes(k)) nlocmin=nlocmin+1
  enddo
  if (nlocmin.gt.0) print*,'nlocmin after',nlocmin

  deallocate (basins,outlets,conn_basins,conn_nodes,conn_weights,active_nodes,tree,sills)

  call loc_min_3_find_donors (rec,ndon,donor,nn)

  call loc_min_3_find_stack (rec,ndon,donor,stack,nn)

  return
end subroutine LocalMinima

!----------------------

subroutine compute_basins (stack,rec,basins,outlets,n,nbasins)

  !Input:

  !stack: parse order from lower to upper nodes
  !receivers: array of recievers

  !Outputs:
  !basins: id of the basin for each node
  !outlets: minimal node of each basin

  !return nbasins: the number of basins

  implicit none

  integer stack(n),rec(n),basins(n),outlets(n)
  integer n,nbasins
  integer ibasins,i,istack,irec

  ibasins=0

  do i=1,n
    istack=stack(i)
    irec=rec(istack)
    if (irec.eq.istack) then
      ibasins=ibasins+1
      outlets(ibasins)=istack
    endif
    basins(istack)=ibasins
    nbasins=ibasins
  enddo

  return
end subroutine compute_basins

!----------------------

subroutine connect_basins (nbasins,basins,outlets,rec,stack,active_nodes,h,nx,ny, &
  basin0,nconn,conn_basins,conn_nodes,conn_weights,nconn_max)
  !Input:
  !nbasins: the number of basins
  !basins: id of the basin for each node
  !outlets: minimal node of each basin
  !receivers: array of recievers
  !stack: parse order from lower to upper nodes
  !active_nodes: nodes not considered as glabal minima (not base nodes)
  !elevation: elevation of each node
  !nx, ny: size of the dem

  !Output:
  !basin0: id of the root basin (the one representing the sea)
  !nconn: number of connections between basins
  !conn_basins: array of pairs of basins (b0, b1) that share a pass
  !conn_nodes: array of pairs (p0, p1) : nodes of the passes for each pair of basin
  !conn_weights : height of the passes (max of elevation[p0], elevation[p1])

  !Connect adjacent basins together through their lowest pass.

  !Creates an (undirected) graph of basins and their connections.

  !The following information is stored for each edge of the graph:

  !- the two grid nodes that togheter form the pass of lowest
  !  elevation found between the two basins ;
  !- a weight that corresponds to the elevation of the pass, i.e.,
  !  the highest elevation among the two nodes forming the pass.

  !Notes
  !-----
  !Connections between open basins are handled differently:

  !Instead of finding connections between adjacent basins,
  !virtual connections are added between one given basin
  !and all other basins.
  !This may save a lot of uneccessary computation, while it
  !ensures a connected graph (i.e., every node has at least
  !an edge), as required for applying minimum spanning tree
  !optimization.

  implicit none

  integer nbasins,basins(nx*ny),outlets(nx*ny),rec(nx*ny),stack(nx*ny)
  integer conn_basins(nconn_max,2),conn_nodes(nconn_max,2)
  double precision conn_weights(nconn_max)
  logical active_nodes(nx*ny)
  double precision h(nx*ny)
  integer nx,ny,n,basin0,nconn,nconn_max

  integer i,j,istack,iused,irec,iistack,iiused,ii(8),jj(8),ki,kj,k
  integer ineighbor,ineighbor_basin,ineighbor_outlet
  integer iconn,ibasin,conn_pos_used_size,conn_idx
  integer, dimension(:), allocatable :: conn_pos,conn_pos_used
  logical active
  double precision weight

  ! theory of planar graph -> max nb. of connections known

  iconn=1
  n=nx*ny

  basin0=-1
  ibasin=1

  ii=(/-1,1,0,0,-1,1,-1,1/)
  jj=(/0,0,-1,1,-1,-1,1,1/)

  allocate (conn_pos(nbasins),conn_pos_used(nbasins))

  conn_pos=-1
  conn_pos_used_size=0

  active=.False.

  ! king (D8) neighbor lookup

  do iistack=1,n
    istack=stack(iistack)
    irec=rec(istack)
    ! new basin
    if (irec.eq.istack) then
      ibasin=basins(istack)
      active=active_nodes(istack)
      do iiused=1,conn_pos_used_size
        iused=conn_pos_used(iiused)
        conn_pos(iused)=-1
      enddo
      conn_pos_used_size=0
      if (.not.active) then
        !      print*,ibasin,basin0
        if (basin0.eq.-1) then
          basin0=ibasin
        else
          conn_basins(iconn,1)=basin0
          conn_basins(iconn,2)=ibasin
          conn_nodes(iconn,1)=-1
          conn_nodes(iconn,2)=-1
          conn_weights(iconn)=-1.d10
          iconn=iconn+1
          !        print*,iconn
        endif
      endif
    endif

    if (active) then

      i=mod(istack-1,nx)+1
      j=(istack-1)/nx+1

      do k=1,8
        ki=i+ii(k)
        kj=j+jj(k)

        if (ki.lt.1.or.ki.gt.nx.or.kj.lt.1.or.kj.gt.ny) goto 111

        ineighbor = (kj-1)*nx+ki
        ineighbor_basin = basins(ineighbor)
        ineighbor_outlet = outlets(ineighbor_basin)

        ! skip same basin or already connected adjacent basin
        ! don't skip adjacent basin if it's an open basin
        if (ibasin.ge.ineighbor_basin.and.active_nodes(ineighbor_outlet)) goto 111

        weight = max(h(istack), h(ineighbor))
        conn_idx = conn_pos(ineighbor_basin)

        ! add new connection
        if (conn_idx .eq. -1) then
          conn_basins(iconn,1) = ibasin
          conn_basins(iconn,2) = ineighbor_basin
          conn_nodes(iconn,1) = istack
          conn_nodes(iconn,2) = ineighbor
          conn_weights(iconn) = weight
          conn_pos(ineighbor_basin) = iconn
          iconn = iconn+1
          conn_pos_used_size = conn_pos_used_size+1
          conn_pos_used(conn_pos_used_size) = ineighbor_basin
          !update existing connection
        elseif (weight .lt. conn_weights(conn_idx)) then
          conn_nodes(conn_idx,1) = istack
          conn_nodes(conn_idx,2) = ineighbor
          conn_weights(conn_idx) = weight
        endif

        111   continue
      enddo

    endif

  enddo

  nconn = iconn - 1

  return
end subroutine connect_basins

!----------------------

subroutine mst_kruskal(conn_weights, conn_basins, nbasins, nconn, mstree, mstree_size)

  !kruskal algorithm to compute a Minimal Spanning Tree

  !Input:
  !conn_basins: array of pairs of basins (b0, b1) that share a pass
  !conn_weights : height of the passes (max of elevation[p0], elevation[p1])
  !nbasins: number of basins

  !Output:
  !mstree: id of the connections in the minimal spanning tree

  implicit none

  double precision conn_weights(nconn)
  integer conn_basins(nconn,2)
  integer mstree(nbasins-1)
  integer nbasins,nconn

  integer, dimension(:), allocatable :: sort_id,parent,rank
  integer mstree_size,eid,eeid,f0,f1,b0,b1

  allocate (sort_id(nconn))
  mstree_size = 0

  ! sort edges
  call loc_min_3_indexx (nconn,conn_weights,sort_id)
  !print*,'weights',conn_weights(sort_id)

  allocate (parent(nbasins),rank(nbasins))
  call UnionFindInit (parent,rank,nbasins)

  do eeid=1,nconn
    eid = sort_id(eeid)

    b0 = conn_basins(eid, 1)
    b1 = conn_basins(eid, 2)

    call Unionfind (b0,parent,nbasins,f0)
    call Unionfind (b1,parent,nbasins,f1)

    if (f0 .ne. f1) then
      mstree_size = mstree_size+1
      mstree(mstree_size) = eid
      call DoUnion (b0,b1,parent,rank,nbasins)
    endif

  enddo

  deallocate (parent,rank)

  return
end subroutine mst_kruskal

!----------------------

subroutine order_tree(edges_n1, edges_n2, edges_p1, edges_p2, nconn, num_nodes, root, &
  tree, ntree, sills, basin_stack, continuous_flow)

  !Swap node order in each edge to follow the inverse of the flow
  !Input:
  !edges_n1, edges_n2:  pairs of basins (b0, b1) that share a pass
  !edges_p1, edges_p2: nodes of the passes for each pair of basin
  !num_nodes: number of basins
  !root: root basin (the sea)
  !tree: tree that need ordering

  !Output:
  !None: edges_n1, edges_n2, edges_p1, edges_p2 are updated in place

  implicit none

  integer num_nodes,root,ntree,nconn
  integer edges_n1(nconn),edges_n2(nconn)
  integer edges_p1(nconn),edges_p2(nconn)
  integer tree(ntree), sills(num_nodes, 2)
  integer basin_stack(num_nodes)
  logical continuous_flow

  integer, dimension(:), allocatable :: nodes_connects_size,nodes_connects_ptr,nodes_adjacency
  integer, dimension(:,:), allocatable :: stack

  integer i,ii,n1,n2,stack_size,nodes_adjacency_size,node,parent,edge_id,edge_swap, basin_stack_size

  !nodes connections
  allocate (nodes_connects_size(num_nodes),nodes_connects_ptr(num_nodes))
  nodes_connects_size=0

  !print*,'num_nodes,ntree,root',num_nodes,ntree,root

  !parse the edges to compute the number of edges per node
  do ii = 1,ntree
    i=tree(ii)
    nodes_connects_size(edges_n1(i)) = nodes_connects_size(edges_n1(i))+1
    nodes_connects_size(edges_n2(i)) = nodes_connects_size(edges_n2(i))+1
  enddo
  !write(*,'(a,1024i4)')'nodes_connect_size',nodes_connects_size(1:num_nodes)

  !compute the id of first edge in adjacency table
  nodes_connects_ptr(1) = 1
  do i=2,num_nodes
    nodes_connects_ptr(i) = nodes_connects_ptr(i-1) + nodes_connects_size(i-1)
    nodes_connects_size(i-1) = 0
  enddo
  !write(*,'(a,1024i4)')'nodes_connect_ptr',nodes_connects_ptr(1:num_nodes)
  !write(*,'(a,1024i4)')'nodes_connect_size',nodes_connects_size(1:num_nodes)

  !create the adjacency table
  nodes_adjacency_size = nodes_connects_ptr(num_nodes) +  nodes_connects_size(num_nodes)-1
  nodes_connects_size(num_nodes) = 0
  !print*,'node_adjacency_size',nodes_adjacency_size
  allocate (nodes_adjacency(nodes_adjacency_size))
  nodes_adjacency=0

  !parse the edges to update the adjacency
  do ii = 1,ntree
    i=tree(ii)
    n1 = edges_n1(i)
    n2 = edges_n2(i)
    nodes_adjacency(nodes_connects_ptr(n1) + nodes_connects_size(n1)) = i
    nodes_adjacency(nodes_connects_ptr(n2) + nodes_connects_size(n2)) = i
    nodes_connects_size(n1) = nodes_connects_size(n1)+1
    nodes_connects_size(n2) = nodes_connects_size(n2)+1
  enddo
  !write(*,'(a,1024i4)')'nodes_connect_ptr',nodes_connects_ptr(1:num_nodes)
  !write(*,'(a,1024i4)')'nodes_connect_size',nodes_connects_size(1:num_nodes)

  !depth-first parse of the tree, starting from root
  !stack of node, parent
  allocate(stack(num_nodes,2))
  stack_size = 1
  stack(1,1) = root
  stack(1,2) = root

  ! sea has no sill
  sills(root,1) = -1
  sills(root,2) = -1

  !print*,'nodes_adjacency',nodes_adjacency
  !print*,'node_ptr',nodes_connects_ptr
  !print*,'nodes_size',nodes_connects_size

  if (continuous_flow) then
    sills(root, 1) = -1
    sills(root, 2) = -1
    basin_stack(1) = root
    basin_stack_size = 1
  endif

  do while (stack_size.gt.0)
    !get parsed node
    node = stack(stack_size, 1)
    parent = stack(stack_size, 2)
    stack_size = stack_size-1
    !  print*,'node',node
    !print*,'range',nodes_connects_ptr(node), nodes_connects_size(node)

    !for each edge of the graph
    !print*,'node',node
    do i=nodes_connects_ptr(node), nodes_connects_ptr(node) + nodes_connects_size(node)-1
      !    print*,i
      edge_id = nodes_adjacency(i)
      !the edge comming from the parent node has already been updated.
    !in this case, the edge is (parent, node)

    if (edges_n1(edge_id).eq.parent.and.node.ne.parent) then

      if (continuous_flow) then
        basin_stack_size = basin_stack_size+1
        basin_stack(basin_stack_size) = node

        sills(node, 1) = edges_p1(edge_id)
        sills(node, 2) = edges_p2(edge_id)
      endif

    else

      !check if the node is not in n1
      if (node.ne.edges_n1(edge_id)) then
        !swap n1 and n2
        edge_swap = edges_n1(edge_id)
        edges_n1(edge_id) = edges_n2(edge_id)
        edges_n2(edge_id) = edge_swap
        !swap p1 and p2
        edge_swap = edges_p1(edge_id)
        edges_p1(edge_id) = edges_p2(edge_id)
        edges_p2(edge_id) = edge_swap
      endif
      !add the opposite node to the stack
      stack_size = stack_size+1
      stack(stack_size,1) =  edges_n2(edge_id)
      stack(stack_size,2) = node

    endif

  enddo

enddo

return
end subroutine order_tree

!----------------------

subroutine correct_receivers(receivers,dist2receivers,outlets,conn_basins,conn_nodes, &
  tree,elevation,nx,ny,dx,dy,nbasins,nconn,ntree)

  !Correct receivers: correct the receivers according the the tree order

  !Input:
  !receivers : array of receviers
  !dist2receivers : distance between bode and its receiver
  !outlets: outlets: minimal node of each basin
  !conn_basins: array of pairs of basins (b0, b1) that share a pass
  !conn_nodes: array of pairs (p0, p1) : nodes of the passes for each pair of basin
  !tree: id of the connections in the minimal spanning tree
  !elevation: elevation of each node
  !nx size of the demin x direction
  !dx, dy: size of a cell

  !Output: none: receivers and dist2receivers are updated in place

  implicit none

  integer nx,ny,nbasins,ntree,nconn
  integer receivers(nx*ny),outlets(nbasins),tree(ntree)
  integer conn_basins(nconn,2),conn_nodes(nconn,2)
  double precision dist2receivers(nx*ny),elevation(nx*ny)
  double precision dx,dy

  integer i,ii,node_from,node_to,outlet_from

  do ii=1,ntree
    i=tree(ii)
    ! tree order: inverse of water flow
    node_to = conn_nodes(i, 1)
    node_from = conn_nodes(i, 2)
    if (node_from.eq.-1) goto 111

    outlet_from = outlets(conn_basins(i, 2))

    ! -> no river erosion on sinks
    dist2receivers(outlet_from) = 1.d10

    if (elevation(node_from).lt.elevation(node_to)) then
      receivers(outlet_from) = node_to
    else
      !receivers[outlet_from] = node_to
      !continue
      receivers(outlet_from) = node_from
      receivers(node_from) = node_to
      !distance based on previous king (4D) neighbor lookup
      if (mod(node_from,nx).eq.mod(node_to,nx)) then
        dist2receivers(node_from) = dx
      else
        dist2receivers(node_from) = dy
      endif

    endif

    111 continue
  enddo

  return
end subroutine correct_receivers

!----------------------

subroutine correct_receivers_v2(receivers,dist2receivers,outlets,conn_basins, &
  conn_nodes,tree,elevation,nx,ny,dx,dy,nbasins,nconn,ntree)

  !Correct receivers: correct the receivers according the the tree order

  !Input:
  !receivers : array of receviers
  !dist2receivers : distance between bode and its receiver
  !outlets: outlets: minimal node of each basin
  !conn_basins: array of pairs of basins (b0, b1) that share a pass
  !conn_nodes: array of pairs (p0, p1) : nodes of the passes for each pair of basin
  !tree: id of the connections in the minimal spanning tree
  !elevation: elevation of each node
  !nx size of the demin x direction
  !dx, dy: size of a cell

  !Output: none: receivers and dist2receivers are updated in place

  implicit none

  integer nx,ny,nbasins,ntree,nconn
  integer receivers(nx*ny),outlets(nbasins),tree(ntree)
  integer conn_basins(nconn,2),conn_nodes(nconn,2)
  integer next_node,cur_node,rcv_next_node
  double precision dist2receivers(nx*ny),elevation(nx*ny)
  double precision dx,dy,ddx,ddy,previous_dist,tmp

  integer i,ii,node_from,node_to,outlet_from

  do ii=1,ntree
    i=tree(ii)
    ! tree order: inverse of water flow
    node_to=conn_nodes(i,1)
    node_from=conn_nodes(i,2)
    if (node_from.eq.-1) goto 111

    outlet_from=outlets(conn_basins(i,2))

    next_node=receivers(node_from)
    cur_node=node_from
    previous_dist=dist2receivers(node_from)

    receivers(node_from)=node_to
    ddx=dx*(mod(node_from-1,nx)-mod(node_to-1,nx))
    ddy=dy*((node_from-1)/nx-(node_to-1)/nx)

    dist2receivers(node_from)=sqrt(ddx*ddx+ddy*ddy)

    do while (cur_node.ne.outlet_from)
      rcv_next_node=receivers(next_node)
      receivers(next_node)=cur_node
      tmp=previous_dist
      previous_dist=dist2receivers(next_node)
      dist2receivers(next_node)=tmp
      cur_node=next_node
      next_node=rcv_next_node
    end do

    111 continue
  enddo

  return
end subroutine correct_receivers_v2


!----------------------

subroutine update_fake_topography ( sills, basin_stack, basins, elevation, nx, ny, dx, dy, nbasins)

  !Input:

  ! sills: pairs of possible sills for each basin
  ! basin_stack: order of basin parsing from sea to top
  ! basins: basin id for each node
  ! nx, ny dx, dy: size and scale of dem


  implicit none

  integer sills(nbasins, 2), basin_stack(nbasins)
  integer basins(nx*ny)
  double precision elevation(nx*ny)
  integer nx,ny,nn,nbasins,i_sill, parse_begin, parse_end, i_node
  integer ii(8),jj(8), i_b, b, nb_dir_i, i_nb, j_nb, i_neighbor
  double precision dx, dy, slope, dist_x, dist_y, new_height

  integer, dimension(:), allocatable :: parse, parsed

  nn = nx*ny

  slope = 1.d-10

  ! parse queue
  ! sills can be parsed twice, so allocate nn*2 for safety
  allocate (parse(2*nn))

  !state of the nodes: "basins" means belong to basin and is not parsed ; -1 means parsed
  allocate (parsed(nn))
  parsed = basins

  !a breadth first parse will be performed, a queue structure is then necessary. "end" means last+1
  parse_begin = 1
  parse_end = 1

  !d8 neighborhood
  ii = (/-1,1,0,0,-1,1,-1,1/)
  jj = (/0,0,-1,1,-1,-1,1,1/)

  !parse basins in flow order
  do i_b = 1, nbasins
    b = basin_stack(i_b)

    !boudary basins have sills -1.
    if (sills(b, 1) .eq. -1) goto 111

    !find sill
    if (elevation(sills(b, 1)) .ge. elevation(sills(b, 2))) then
      i_sill = sills(b, 1)
    else
      i_sill = sills(b, 2)
    endif

    ! start parsing at sill
    parse(parse_begin) = i_sill
    parse_end = parse_end+1

    parsed(i_sill) = -1

    do while (parse_begin.ne.parse_end)

      i_node = parse(parse_begin)
      parse_begin = parse_begin+1

      ! parse neighbors
      do nb_dir_i=1,8

        i_nb = mod(i_node-1,nx)+1 + ii(nb_dir_i)
        j_nb = (i_node-1)/nx+1 + jj(nb_dir_i)

        if (i_nb.ge.1 .and. i_nb.le.nx .and. j_nb.ge.1 .and. j_nb.le.ny) then
          i_neighbor = (j_nb-1) * nx + i_nb
          if (parsed(i_neighbor) .eq. b) then
            parsed(i_neighbor) = -1
            dist_x = ii(nb_dir_i) * dx
            dist_y = jj(nb_dir_i) * dy
            new_height = elevation(i_node) + slope * sqrt(dist_x*dist_x + dist_y*dist_y)
            if (elevation(i_neighbor) .le. new_height) then
              elevation(i_neighbor) = new_height
              parse(parse_end) = i_neighbor
              parse_end = parse_end+1
            endif
          endif
        endif

      enddo

    enddo

    111 continue
  enddo

  return
end subroutine update_fake_topography

!----------------------

subroutine UnionFindInit (parent,rank,n)

  implicit none

  integer parent(n),rank(n)
  integer n,i

  do i=1,n
    parent(i)=i
    rank(i)=0
  enddo

  return
end subroutine UnionFindInit

!----------------------

subroutine DoUnion (x,y,parent,rank,n)

  implicit none

  integer parent(n),rank(n)
  integer n,x,y
  integer xroot,yroot

  call UnionFind(x,parent,n,xroot)
  call UnionFind(y,parent,n,yroot)

  if (xroot.ne.yroot) then

    if (rank(xroot).lt.rank(yroot)) then
      parent(xroot)=yroot

    else
      parent(yroot)=xroot

      if (rank(xroot).eq.rank(yroot)) then
        rank(xroot)=rank(xroot)+1
      endif

    endif

  endif

  return
end subroutine DoUnion

!----------------------

subroutine UnionFind (x,parent,n,found)

  implicit none

  integer x,xp,xc,n,found
  integer parent(n)

  xp=x
  xc=-1000
  do while (xc.ne.xp)
    xc=xp
    xp=parent(xc)
  enddo
  parent(x)=xc
  found=xc

  return
end subroutine UnionFind

!----------------------

subroutine loc_min_3_find_receivers (h,rec,length,bc,nx,ny,dx,dy)

  implicit none

  double precision h(nx*ny),length(nx*ny),dx,dy
  integer rec(nx*ny)
  logical bc(nx*ny)
  integer nx,ny

  integer i,j,ij,ii,jj,iii,jjj,ijk
  double precision smax,l,slope

  do j=1,ny
    do i=1,nx
      ij=i+(j-1)*nx
      smax=tiny(smax)
      rec(ij)=ij
      if (bc(ij)) goto 111
      do jj=-1,1
        do ii=-1,1
          iii=i+ii
          jjj=j+jj
          if (iii.ge.1.and.jjj.ge.1.and.iii.le.nx.and.jjj.le.ny) then
            ijk=iii+(jjj-1)*nx
            if (ijk.ne.ij) then
              l=sqrt((dx*ii)**2+(dy*jj)**2)
              slope=(h(ij)-h(ijk))/l
              if (slope.gt.smax) then
                smax=slope
                rec(ij)=ijk
                length(ij)=l
              endif
            endif
          endif
        enddo
      enddo
      111 continue
    enddo
  enddo

  return
end subroutine loc_min_3_find_receivers

!----------------------

subroutine loc_min_3_find_donors (rec,ndon,donor,nn)

  implicit none

  integer nn
  integer rec(nn),ndon(nn),donor(8,nn)

  integer ij,ijk

  ndon=0
  do ij=1,nn
    if (rec(ij).ne.ij) then
      ijk=rec(ij)
      ndon(ijk)=ndon(ijk)+1
      donor(ndon(ijk),ijk)=ij
    endif
  enddo

  return
end subroutine loc_min_3_find_donors

!----------------------

subroutine loc_min_3_find_stack (rec,ndon,donor,stack,nn)

  implicit none

  integer nn
  integer rec(nn),stack(nn),ndon(nn),donor(8,nn)

  integer nstack,ij

  nstack=0
  do ij=1,nn
    if (rec(ij).eq.ij) then
      nstack=nstack+1
      stack(nstack)=ij
      call find_stack_recursively_locmin (ij,donor,ndon,nn,stack,nstack)
    endif
  enddo

  return
end subroutine loc_min_3_find_stack

!----------------------

recursive subroutine loc_min_3_find_stack_recursively  (ij,donor,ndon,nn,stack,nstack)

implicit none

integer ij,nn,nstack
integer donor(8,nn),ndon(nn),stack(nn)

integer k,ijk

do k=1,ndon(ij)
  ijk=donor(k,ij)
  nstack=nstack+1
  stack(nstack)=ijk
  call loc_min_3_find_stack_recursively (ijk,donor,ndon,nn,stack,nstack)
enddo

return
end subroutine loc_min_3_find_stack_recursively

!-----------------------

subroutine loc_min_3_indexx(n,arr,indx)

  implicit none

  integer n,indx(n),M,NSTACK
  double precision arr(n)
  parameter (M=7,NSTACK=50)
  integer i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
  double precision a

  do j=1,n
    indx(j)=j
  enddo

  jstack=0
  l=1
  ir=n

  1 if(ir-l.lt.M)then

    do j=l+1,ir
      indxt=indx(j)
      a=arr(indxt)
      do i=j-1,1,-1
        if(arr(indx(i)).le.a)goto 2
        indx(i+1)=indx(i)
      enddo
      i=0
      2   indx(i+1)=indxt
    enddo

    if(jstack.eq.0) return
    ir=istack(jstack)
    l=istack(jstack-1)
    jstack=jstack-2

  else

    k=(l+ir)/2
    itemp=indx(k)
    indx(k)=indx(l+1)
    indx(l+1)=itemp

    if(arr(indx(l+1)).gt.arr(indx(ir)))then
      itemp=indx(l+1)
      indx(l+1)=indx(ir)
      indx(ir)=itemp
    endif

    if(arr(indx(l)).gt.arr(indx(ir)))then
      itemp=indx(l)
      indx(l)=indx(ir)
      indx(ir)=itemp
    endif

    if(arr(indx(l+1)).gt.arr(indx(l)))then
      itemp=indx(l+1)
      indx(l+1)=indx(l)
      indx(l)=itemp
    endif

    i=l+1
    j=ir
    indxt=indx(l)
    a=arr(indxt)
    3 continue
    i=i+1
    if(arr(indx(i)).lt.a) goto 3

    4 continue
    j=j-1
    if(arr(indx(j)).gt.a) goto 4

    if(j.lt.i) goto 5
    itemp=indx(i)
    indx(i)=indx(j)
    indx(j)=itemp
    goto 3

    5  indx(l)=indx(j)
    indx(j)=indxt
    jstack=jstack+2

    if(jstack.gt.NSTACK) stop 'NSTACK too small in loc_min_3_indexx'

    if(ir-i+1.ge.j-l)then
      istack(jstack)=ir
      istack(jstack-1)=i
      ir=j-1
    else
      istack(jstack)=j-1
      istack(jstack-1)=l
      l=i
    endif

  endif

  goto 1

end subroutine loc_min_3_indexx

!----------------------

recursive subroutine find_stack_recursively_locmin (ij,donor,ndon,nn,stack,nstack)

implicit none

integer ij,nn,nstack
integer donor(8,nn),ndon(nn),stack(nn)

integer k,ijk

do k=1,ndon(ij)
  ijk=donor(k,ij)
  nstack=nstack+1
  stack(nstack)=ijk
  call find_stack_recursively_locmin (ijk,donor,ndon,nn,stack,nstack)
enddo

return
end subroutine find_stack_recursively_locmin
