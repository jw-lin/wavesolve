import pygmsh
import numpy as np
import matplotlib.pyplot as plt
import gmsh
import math

#region miscellaneous functions   

def get_unique_edges(mesh,mutate=True):
    """ get set of unique edges in mesh. when mutate = True, this function adds a connectivity table for edges to the mesh;
        this table is not included in the base pygmsh/gmsh meshes, but is needed for vectorial FEM. """
    tris = mesh.cells[1].data
    
    unique_edges = list()
    edge_indices = np.zeros((tris.shape[0],3),dtype=np.uint)
    edge_flips = np.ones((tris.shape[0],3))

    i = 0
    for j,tri in enumerate(tris):
        e0 = sorted([tri[0],tri[1]])
        e1 = sorted([tri[1],tri[2]])
        e2 = sorted([tri[2],tri[0]])
        es = [e0,e1,e2]
        es_unsort = [[tri[0],tri[1]],[tri[1],tri[2]],[tri[2],tri[0]]]
        for k,e in enumerate(es):
            if e not in unique_edges:
                unique_edges.append(e)
                edge_indices[j,k] = i
                i += 1
            else:
                edge_indices[j,k] = unique_edges.index(e)
            if e != es_unsort[k]:
                edge_flips[j,k] = -1

    out = np.array(unique_edges)
    if mutate:
        mesh.cells[0].data = out
        mesh.edge_indices = edge_indices
        mesh.edge_flips = edge_flips
        mesh.num_edges = len(out)
        mesh.num_points = mesh.points.shape[0]
    return out,edge_indices

def intersect(a,b,c,d):
    # credit to stackoverflow guy
    a1 = b[1]-a[1]
    b1 = a[0]-b[0]
    c1 = a1*a[0] + b1*a[1]
    a2 = d[1]-c[1]
    b2 = c[0]-d[0]
    c2 = a2*c[0] + b2*c[1]
    det = a1*b2 - a2*b1
    if (det== 0):
        return np.inf, np.inf
    else:
        x = (b2*c1 - b1*c2)/det
        y = (a1*c2 - a2*c1)/det
        return x,y

def pt_edge_dist(p,a,b):
    """ compute distance between point p and line segment with points a and b.
        a and b may also be Nx2 arrays corresponding to N edges.
    """
    if a.ndim == 1:
        ab = b - a
        ap = p - a
        squared_length_ab = np.dot(ab, ab)
        if squared_length_ab == 0:
            return np.linalg.norm(p - a)
        t = np.dot(ap, ab) / squared_length_ab
        closest = a + np.clip(t,0,1) * ab
        distance = np.linalg.norm(p - closest)
        return distance
    else:
        ab = b - a
        ap = p[None,:] - a
        squared_length_ab = np.sum(ab*ab,axis=1)
        t = np.sum(ap*ab,axis=1) / squared_length_ab
        closest = a + np.clip(t,0,1)[:,None] * ab
        distance = np.linalg.norm(p[None,:] - closest,axis=1)
        return distance
    
def plot_mesh(mesh,IOR_dict=None,alpha=0.3,ax=None,plot_points=True):
    """ plot a mesh and associated refractive index distribution
    Args:
    mesh: the mesh to be plotted. if None, we auto-compute a mesh using default values
    IOR_dict: dictionary that assigns each named region in the mesh to a refractive index value
    """
    show=False
    verts=3
    if ax is None:
        fig,ax = plt.subplots(figsize=(5,5))
        ax.set_aspect("equal")
        show=True

    ax.set_aspect('equal')

    points = mesh.points
    els = mesh.cells[1].data
    materials = mesh.cell_sets.keys()

    if IOR_dict is not None:
        IORs = [ior[1] for ior in IOR_dict.items()]
        n,n0 = max(IORs) , min(IORs)

    for material in materials:   
        if IOR_dict is not None:    
            cval = IOR_dict[material]/(n-n0) - n0/(n-n0)
            cm = plt.get_cmap("inferno")
            color = cm(cval)
        else:
            color="None"

        _els = els[tuple(mesh.cell_sets[material])][0,:,0,:]
        for i,_el in enumerate(_els):
            t=plt.Polygon(points[_el[:verts]][:,:2], facecolor=color)
            t_edge=plt.Polygon(points[_el[:verts]][:,:2], lw=0.5,color='0.5',alpha=alpha,fill=False)
            ax.add_patch(t)
            ax.add_patch(t_edge)

    ax.set_xlim(np.min(points[:,0]),np.max(points[:,0]) )
    ax.set_ylim(np.min(points[:,1]),np.max(points[:,1]) )
    if plot_points:
        for point in points:
            ax.plot(point[0],point[1],color='0.5',marker='o',ms=1.5,alpha=alpha)

    if show:
        plt.show()

def boolean_fragment(geom:pygmsh.occ.Geometry,_object,tool):
    """ fragment the tool and the object, and return the fragments in the following order:
        intersection, object_fragment, tool_fragment.
        in some cases one of the later two may be empty
    """
    object_copy = geom.copy(_object)
    tool_copy = geom.copy(tool)
    try:
        intersection = geom.boolean_intersection([object_copy,tool_copy])
    except:
        # no intersection - make first element None to signal
        return [None,_object,tool]

    _object = geom.boolean_difference(_object,intersection,delete_first=True,delete_other=False)
    tool = geom.boolean_difference(tool,intersection,delete_first=True,delete_other=False)
    return intersection+_object+tool

def boolean_difference(geom,_object,_tool):
    if type(_object) == list:
        for o in _object:
            if type(_tool) == list:
                for t in _tool:
                    o = geom.boolean_difference(o,t,delete_other=False,delete_first=True)
            else:
                o = geom.boolean_difference(o,_tool,delete_other=False,delete_first=True)
    else:
        if type(_tool) == list:
            for t in _tool:
                _object = geom.boolean_difference(_object,t,delete_other=False,delete_first=True)
        else:
            _object = geom.boolean_difference(_object,_tool,delete_other=False,delete_first=True)
    
    return _object

def dist(p1,p2):
    return np.sqrt(np.sum(np.power(p1-p2,2)))

def ellipse_dist(semi_major, semi_minor, c, p, iters=3):
    """ compute signed distance to axis-aligned ellipse boundary """  

    _p = [p[0]-c[0],p[1]-c[1]]

    px = abs(_p[0])
    py = abs(_p[1])

    tx = 0.707
    ty = 0.707

    a = semi_major
    b = semi_minor

    inside = _p[0]**2/semi_major**2 + _p[1]**2/semi_minor**2 <= 1

    for x in range(0,iters):
        x = a * tx
        y = b * ty

        ex = (a*a - b*b) * tx**3 / a
        ey = (b*b - a*a) * ty**3 / b

        rx = x - ex
        ry = y - ey

        qx = px - ex
        qy = py - ey

        r = math.hypot(ry, rx)
        q = math.hypot(qy, qx)
        
        if q == 0:
            tx = 1
            ty = 1
        else:
            tx = min(1, max(0, (qx * r / q + ex) / a))
            ty = min(1, max(0, (qy * r / q + ey) / b))
        
        t = math.hypot(ty, tx)
        tx /= t 
        ty /= t 

    isect = [math.copysign(a * tx,_p[0]), math.copysign(b * ty,_p[1])]
    sgn = 1 if not inside else -1
    dist = math.sqrt((_p[0]-isect[0])**2 + (_p[1]-isect[1])**2)*sgn
    return dist

#endregion    

#region Prim2D
class Prim2D:
    """ a Prim2D (2D primitive) is an an array of N (x,y) points, shape (N,2), that denote a closed curve (so, a polygon). 
        inside the closed curve, the primitive has refractive index n. 
    """
    def __init__(self,n,label,points=[]):
        """
        ARGS:
            n (float): refractive index
            label (str): a bookkeeping label to attach to this object (e.g. "core")
            points: a list or array of [x,y] points corresponding to the shape's boundary
        """
        self.points = points
        self.label = label
        self.n = n
        self.res = len(points)
        self.mesh_size = None # set to a numeric value to force a triangle size within the closed region
        self.skip_refinement = False
    
    def __str__(self):
        return type(self).__name__+"('"+self.label+"')"

    def make_poly(self,geom):
        """ convert self.points into a Gmsh polygon

        ARGS:
            geom: pygmsh geometry kernel
        """
        if hasattr(self.points[0][0],'__len__'):
            ps = [geom.add_polygon(p) for p in self.points]
            poly = geom.boolean_union(ps)[0]
        else:
            poly = geom.add_polygon(self.points)
        return poly

    def make_points(self,args):
        """ make an Nx2 array of points for marking the primitive boundary,
            according to some args. CCW ordering! this is a placeholder; 
            subclasses should implement this function.
        """
        raise NotImplementedError         

    def inside(self,x,y):
        """
        ARGS:
            x (float): x coordinate of a point
            y (float): y coordinate of a point
        
        RETURNS:
            (bool) : whether or not [x,y] is inside the primitive boundary
        """
        raise NotImplementedError
    
    def sgn_inside(self,x,y):
        if self.inside(x,y):
            return -1
        return 1

    def boundary_dist(self,x,y):
        """ compute the signed distance of the point (x,y) to the primitive boundary.
            this default implementation treats the primitive as a polygon and computes
            the point-edge distance wrt the every edge, then takes the minimum. 
            it relies on inside() to convert distance to a signed distance.

        ARGS:
            x (float): x coordinate of a point
            y (float): y coordinate of a point
        
        RETURNS:
            (float) : signed distance to the boundary (negative inside, positive outside)
        """
        if hasattr(self,"boundary_pts"):
            points = self.boundary_pts
        else:
            points = self.points
        
        pts1 = points
        pts2 = np.roll(points,1,0)
        dists = pt_edge_dist(np.array([x,y]),pts1,pts2)
        #dists = np.array([pt_edge_dist(np.array([x,y]),pts1[i],pts2[i]) for i in range(points.shape[0])])  #
        return np.min(dists)*self.sgn_inside(x,y)
    
    def dist_map(self,bounds,res=50):
        """ plot a signed distance map. useful for checking boundary_dist()
        
        ARGS:
            bounds (array): [xmin,xmax,ymin,ymax] extent of the plot
            res (int): 1D resolution of the plot
        """
        xa = np.linspace(bounds[0],bounds[1],res)
        ya = np.linspace(bounds[2],bounds[3],res)
        dists = np.empty((res,res))
        fig,ax = plt.subplots(1,1)
        ax.set_aspect('equal')
        for i in range(res):
            for j in range(res):
                dists[i,j] = self.boundary_dist(xa[i],ya[j])
        im = ax.imshow(dists.T,origin="lower",extent=bounds)
        self.plot_boundary(ax)
        plt.colorbar(im,ax=ax)
        plt.show()

    def get_nearest_bp_idx(self,point):
        """ get the nearest boundary point index for a given point """
        dists = [np.abs(dist(point,p)) for p in self.points]
        return np.argmin(dists)

    def plot_boundary(self,ax=None,color="red",alpha=1,lw=1):
        """ plot the boundary of the primitive 

        ARGS:
            ax (matplotlib.axes): matplotlib axis to plot on; if None, one will be made.
            color (str): linecolor of the boundary
            alpha (float): line opacity
            lw (float): line width
        """
        if hasattr(self,"boundary_pts"):
            points = self.boundary_pts
        else:
            points = self.points
        show = False
        if ax is None:
            show = True
            fig,ax = plt.subplots(1,1)
            ax.set_aspect("equal")
        points_arr =  np.zeros((points.shape[0]+1,points.shape[1]))
        points_arr[:-1] = points
        points_arr[-1] = points_arr[0]
        ax.plot(points_arr.T[0],points_arr.T[1],color=color,alpha=alpha,lw=lw)
        if show:
            plt.show() 

class Circle(Prim2D):
    """ a (discretized) Circle primitive, defined by radius, center, and number of sides """

    def make_points(self,radius,res,center=(0,0)):
        """
        ARGS:
            radius (float): circle radius
            res (int): number of line segments with which to discretize the circle
            center: (x,y) coordinates of the circle center
        """
        thetas = np.linspace(0,2*np.pi,res,endpoint=False)
        points = []
        for t in thetas:
            points.append((radius*np.cos(t)+center[0],radius*np.sin(t)+center[1]))        
        points = np.array(points)

        self.radius = radius # save params for later comp
        self.center = center # 
        self.points = points
        return points
    
    def inside(self, x, y):
        return self.boundary_dist(x,y) < 0

    def boundary_dist(self, x, y):
        return np.sqrt(np.power(x-self.center[0],2)+np.power(y-self.center[1],2)) - self.radius

class Rectangle(Prim2D):
    """ rectangle primitive, defined by corner pounts. """

    def make_points(self,xmin,xmax,ymin,ymax):
        """
        ARGS:
            xmin (float): minimum x coordinate of the rectangle
            xmax (float): maximum x coordinate
            ymin (float): minimum y coordinate
            ymax (float): maximum y coordinate
        """
        points = np.array([[xmin,ymin],[xmax,ymin],[xmax,ymax],[xmin,ymax]])
        self.bounds = [xmin,xmax,ymin,ymax]
        self.points = points
        self.center = [(xmin+xmax)/2,(ymin+ymax)/2]
        return points

    def inside(self,x,y):
        if (self.bounds[0]<=x<=self.bounds[1]) and (self.bounds[2]<=y<=self.bounds[3]):
            return True
        return False

    def boundary_dist(self,x,y):
        # i got this from gemini ...
        p = np.array([x,y])
        c = np.array(self.center)
        sz = np.array([self.bounds[1]-self.bounds[0],self.bounds[3]-self.bounds[2]])/2.
        d = np.abs(p - c) - sz
        outside_dist = np.linalg.norm(np.maximum(d, 0.0))
        inside_dist = np.maximum(d[0], d[1]) if np.any(d < 0) else 0.0
        return outside_dist + np.minimum(inside_dist, 0.0)

class Ellipse(Prim2D):
    """axis-aligned ellipse. a is the semi-axis along x, b is the semi-axis along y."""
    def make_points(self,a,b,res,center=(0,0)):
        """
        ARGS:
            a (float): semimajor axis
            b (float): semiminor axis
            res (int): number of line segments with which to discretize the ellipse
            center: (x,y) coordinates of the circle center
        """
        thetas = np.linspace(0,2*np.pi,res,endpoint=False)
        points = []
        for t in thetas:
            points.append((a*np.cos(t)+center[0],b*np.sin(t)+center[1]))        
        points = np.array(points)

        self.a = a
        self.b = b
        self.center = center
        self.points = points
        return points

    def inside(self, x, y):
        return self.boundary_dist(x,y) < 0

    def boundary_dist(self, x, y):
        return ellipse_dist(self.a,self.b,self.center,[x,y])

class Prim2DUnion(Prim2D):
    """ a union of Prim2Ds """
    def __init__(self,ps:list[Prim2D],label):
        """
        ARGS:
            ps (list[Prim2D]): a list of Prim2D objects to boolean union
            label (str): a bookkeeping label to attach to this object (e.g. "core")
        """
        ns = [p.n for p in ps]
        assert np.all(np.array(ns)==ns[0]),"primitives must have the same refractive index"
        centers = np.array([p.center for p in ps])
        self.center = np.mean(centers,axis=0)
        points = [p.points for p in ps]
        super().__init__(ps[0].n,label,points)
        self.ps = ps
        self.boundary_pts = self.get_boundary_points()

    def make_points(self,args):
        out = []
        for i,p in enumerate(self.ps):
            points = p.make_points(args[i])
            out.append(points)
        return out
    
    def inside(self,x,y):
        for p in self.ps:
            if p.inside(x,y):
                return True
        return False

    def make_poly(self,geom):
        if hasattr(self.points[0][0],'__len__'):
            ps = [geom.add_polygon(p) for p in self.points]
            polys = geom.boolean_union(ps)
            poly = polys
        else:
            poly = geom.add_polygon(self.points)
        return poly
    
    def get_boundary_points(self):
        # this was a nightmare to write
        ps = self.ps
        max_dists = [ np.max(np.linalg.norm(p.points,axis=1)) for p in ps ]
        pts_idx = np.argmax(max_dists)
        pt_idx = np.argmax(np.linalg.norm(ps[pts_idx].points,axis=1))
        eps=1e-10
        boundary_pts = []
        test_point = ps[pts_idx].points[pt_idx]
        Np = len(ps[pts_idx].points)
        counter = 0
        while (len(boundary_pts) == 0 or dist(test_point,boundary_pts[0]) > eps):
            dists = []
            for j,p in enumerate(ps):
                if j == pts_idx:
                    dists.append(np.inf)
                    continue
                d = p.boundary_dist(test_point[0],test_point[1])
                dists.append(d)
            dists = np.array(dists)
            imin = np.argmin(dists)
            dmin = dists[imin]    
            if dmin < 0:
                #we've crossed something
                p_last = ps[pts_idx]
                
                # goto another prim
                pts_idx = imin

                Np = len(ps[pts_idx].points)

                a = boundary_pts[-1]
                b = test_point

                pt_idx = ps[pts_idx].get_nearest_bp_idx(test_point)
                closest_on_new_prim = ps[pts_idx].points[pt_idx]

                if p_last.boundary_dist(closest_on_new_prim[0],closest_on_new_prim[1])<0:
                    c = ps[pts_idx].points[pt_idx]
                    d =  ps[pts_idx].points[(pt_idx+1)%Np]
                    pt_idx = (pt_idx + 1)%Np
                else:
                    c = ps[pts_idx].points[(pt_idx-1)%Np]
                    d =  ps[pts_idx].points[pt_idx]
                isect = intersect(a,b,c,d)
                boundary_pts.append(isect)
                test_point = ps[pts_idx].points[pt_idx]
            else:
                boundary_pts.append(test_point)
                pt_idx = (pt_idx + 1)%Np
                test_point = ps[pts_idx].points[pt_idx]
            counter += 1
        return np.array(boundary_pts)

class Prim2DArray(Prim2D):
    """ an array of non-intersecting Prim2Ds which share the same material.
        this is functionally equivalent to a list of Prim2D objects, except
        that all the Prim2D objects now share a label.
    """
    def __init__(self,ps:list[Prim2D],label):
        """
        ARGS:
            ps (list[Prim2D]): a list of Prim2D objects
            label (str): a bookkeeping label to attach to this object (e.g. "core")
        """
        ns = [p.n for p in ps]
        assert np.all(np.array(ns)==ns[0]),"primitives must have the same refractive index"
        super().__init__(ps[0].n,label,[p.points for p in ps])
        self.ps = ps
        self.label = label
    
    def boundary_dist(self, x, y):
        dists = [p.boundary_dist(x,y) for p in self.ps]
        return np.min(dists)

    def make_poly(self,geom):
        polys = [geom.add_polygon(p) for p in self.points]
        return polys

    def plot_boundary(self, ax=None, color="red", alpha=1, lw=1):
        show = False
        if ax is None:
            show = True
            fig,ax = plt.subplots(1,1)
            ax.set_aspect("equal")
        for p in self.ps:
            p.plot_boundary(ax, color, alpha, lw)
        if show:
            plt.show()

#endregion    

#region Waveguide
        
class Waveguide:
    """ a Waveguide is a collection of prim2Ds, organized into layers. the refractive index 
    of earlier layers is overwritten by later layers.
    """

    #: float: mesh boundary refinement linear distance scaling  
    mesh_dist_scale = 0.25   

    #: float: mesh boundary refinement power scaling
    mesh_dist_power = 1.0

    #: float: minimum allowed mesh size
    min_mesh_size = 0.1
    
    #: float: maximum allowed mesh size
    max_mesh_size = 10.

    def __init__(self,prim2Ds):
        """
        ARGS:
            prim2Ds: a (potentially) nested list of Prim2D objects. Later
                          elements overwrite earlier ones, in terms
                          of refractive index.
        """

        self.prim2Dgroups = prim2Ds # an arrangement of Prim2D objects, stored as a (potentially nested) list. each element is overwritten by the next.
        self.IOR_dict = {}
        
        primsflat = [] # flat array of primitives
        for i,p in enumerate(self.prim2Dgroups):
            if type(p) == list:    
                for _p in p:
                    primsflat.append(_p)
            else:
                primsflat.append(p)  
        
        self.primsflat = primsflat
        self.primsflat[0].skip_refinement = True

    def make_mesh(self,algo=6,order=2,adaptive=True):
        """ Construct a mesh with boundary refinement at material interfaces

        ARGS:
            algo (int): Gmsh algorithm, values 1-7
            order (int): finite element mesh order. 1 -> linear triangle elements, 2 (default) -> quadratic triangle elements
            adaptive (bool): whether to adaptively vary mesh triangle size, default True.
        """

        _scale = self.mesh_dist_scale
        _power = self.mesh_dist_power
        min_mesh_size = self.min_mesh_size
        max_mesh_size = self.max_mesh_size

        with pygmsh.occ.Geometry() as geom:
            gmsh.option.setNumber('General.Terminal', 0)
            # make the polygons
            nested_polygons = []
            for el in self.prim2Dgroups:
                if type(el) != list:
                    #polygons.append(geom.add_polygon(el.prim2D.points))
                    poly = el.make_poly(geom)
                    nested_polygons.append(poly)
                else:
                    els = []
                    nested_els = []
                    for _el in el:
                        poly = _el.make_poly(geom)
                        nested_els.append(poly)
                        if type(poly) == list:
                            els += poly
                        else:
                            els.append(poly)
                    nested_polygons.append(nested_els)

            # diff the polygons
            for i in range(0,len(self.prim2Dgroups)-1):
                polys = nested_polygons[i]
                for j in range(i+1,len(self.prim2Dgroups)):
                    _polys = nested_polygons[j]
                    polys = boolean_difference(geom,polys,_polys)

            # add physical groups
            for i,prim in enumerate(self.prim2Dgroups):
                if type(prim) == list:
                    for j,pprim in enumerate(prim):
                        geom.add_physical(nested_polygons[i][j],pprim.label)
                else:
                    geom.add_physical(nested_polygons[i],prim.label)

            if adaptive:
                # mesh refinement callback
                def callback(dim,tag,x,y,z,lc):
                    return self.compute_mesh_size(x,y,_scale=_scale,_power=_power,min_size=min_mesh_size,max_size=max_mesh_size)
                geom.set_mesh_size_callback(callback)

            geom.env.removeAllDuplicates()
            mesh = geom.generate_mesh(dim=2,order=order,algorithm=algo)
            get_unique_edges(mesh)

            return mesh
        
    def compute_mesh_size(self,x,y,_scale=1.,_power=1.,min_size=None,max_size=None):
        """ compute a target mesh size (triangle side length) at the point (x,y)

        ARGS:
            x (float): x point to compute mesh size at
            y (float): y point to compute mesh size at
            _scale (float): a factor that determines how quickly mesh size should increase away from primitive boundaries. higher = more quickly.
            _power (float): another factor that determines how mesh size increases away from boundaries. default = 1 (linear increase). higher = more quickly.
            min_size (float): the minimum mesh size that the algorithm can choose
            max_size (float): the maximum mesh size that the algorithm can chooose
        """
        prims = self.primsflat
        dists = np.zeros(len(prims)) # compute a distance to each primitive boundary
        for i,p in enumerate(prims): 
            if p.skip_refinement and p.mesh_size is not None:
                dists[i] = 0. # if there is a set mesh size and we dont care about boundary refinement, set dist=0 -> fixed mesh size inside primitive later
            else:
                dists[i] = p.boundary_dist(x,y)
        # compute target mesh sizes
        mesh_sizes = np.zeros(len(prims))
        for i,d in enumerate(dists): 
            p = prims[i]
            if hasattr(p,'boundary_pts'):
                boundary_mesh_size = dist(p.boundary_pts[0],p.boundary_pts[1])
            else:
                boundary_mesh_size = dist(p.points[0],p.points[1])

            if p.mesh_size is not None:
                boundary_mesh_size = min(boundary_mesh_size,p.mesh_size)  

            scaled_size = (1+np.power(np.abs(d)/boundary_mesh_size *_scale ,_power)) * boundary_mesh_size # this goes to boundary_mesh_size as d->0, and increases as d->inf for _power>0
            if d<=0 and p.mesh_size is not None:
                mesh_sizes[i] = min(scaled_size,p.mesh_size)
            else:
                mesh_sizes[i] = scaled_size
        target_size = np.min(mesh_sizes)
        if min_size:
            target_size = max(min_size,target_size)
        if max_size:
            scaled_size = min(max_size,target_size)    
        return target_size

    def assign_IOR(self):
        """ returns a dictionary which maps all material labels in the Waveguide mesh
            to the corresponding refractive index value. """
        for p in self.prim2Dgroups:
            if type(p) == list:
                for _p in p:
                    if _p.label in self.IOR_dict:
                        continue
                    self.IOR_dict[_p.label] = _p.n
            else:
                if p.label in self.IOR_dict:
                    continue
                self.IOR_dict[p.label] = p.n  
        return self.IOR_dict

    def plot_mesh(self,mesh=None,IOR_dict=None,alpha=0.3,ax=None,plot_points=True):
        """ plot a mesh and associated refractive index distribution

        ARGS:
            mesh: the finite element mesh to be plotted. if None, one is made with make_mesh()
            IOR_dict (dict or None): dictionary that assigns each named region in the mesh to a refractive index value. 
                                     if None, one is made automatically.
            alpha (float): opacity of mesh lines
            ax (matplotlib.axes or None): matplotlib axis to plot on; if None, one will be made
            plot_points (bool): plot the mesh triangle points (vertices, and midpoints if order 2)
        """
        if mesh is None:
            mesh = self.make_mesh()
        if IOR_dict is None:
            IOR_dict = self.assign_IOR()

        plot_mesh(mesh,IOR_dict,alpha,ax,plot_points)
    
    @staticmethod
    def plot_boundaries_recursive(prim,ax,color='red',lw=1,alpha=1):
        if type(prim) is list:
            for p in prim:
                Waveguide.plot_boundaries_recursive(p,ax,color,lw,alpha)
        else:
            prim.plot_boundary(ax,color,lw,alpha)

    def plot_boundaries(self,ax=None,color='red',lw=1,alpha=1):
        """ plot the boundaries of all prim2Dgroups

        ARGS:
            ax (matplotlib.axes or None): matplotlib axis to plot on; if None, one will be made
            color (str): plotline color, default "red"
            lw (float): linewidth, default 1
            alpha (float): opacity, default 1
        """
        
        show=False
        if ax is None:
            fig,ax = plt.subplots(figsize=(5,5))
            ax.set_aspect("equal")
            show=True
        else:
            ax.autoscale(False)
        for prim in self.prim2Dgroups:
            self.plot_boundaries_recursive(prim,ax,color,lw,alpha)
        if show:
            plt.show()
    
class CircularFiber(Waveguide):
    """a circular step-index optical fiber"""
    def __init__(self,rcore,rclad,ncore,nclad,core_res,clad_res=None,core_mesh_size=None,clad_mesh_size=None):
        """
        ARGS:
            rcore (float): radius of core
            rclad (float): radius of cladding
            ncore (float): core index
            nclad (float): cladding index
            core_res (int): number of line segments to divide the core boundary into
            clad_res (int): number of line segments to divide the cladding boundary into, default core_res/2
        """
        if clad_res == None:
            clad_res = int(core_res/2)
        core = Circle(ncore,"core")
        core.make_points(rcore,core_res)
        core.mesh_size = core_mesh_size

        cladding = Circle(nclad,"cladding")
        cladding.make_points(rclad,clad_res)
        cladding.mesh_size = clad_mesh_size
        super().__init__([cladding,core])

class EllipticalFiber(Waveguide):
    """an axis-aligned elliptical core step-index fiber"""
    def __init__(self,acore,bcore,rclad,ncore,nclad,core_res,clad_res=None,core_mesh_size=None,clad_mesh_size=None):
        """
        ARGS:
            acore (float): extent of elliptical core along x (the "x" radius)
            bcore (float): extent of elliptical core along y (the "y" radius)
            rclad (float): radius of cladding, assumed circular
            ncore (float): core index
            nclad (float): cladding index
            core_res (int): number of line segments to divide the core boundary into
            clad_res (int): number of line segments to divide the cladding boundary into, default core_res/2
        """
        if clad_res == None:
            clad_res = int(core_res/2)
        core = Ellipse(ncore,"core")
        core.make_points(acore,bcore,core_res)
        core.mesh_size = core_mesh_size
        cladding = Circle(nclad,"cladding")
        cladding.make_points(rclad,clad_res)
        cladding.mesh_size = clad_mesh_size
        super().__init__([cladding,core])

class PhotonicCrystalFiber(Waveguide):
    """ an optical fiber filled with a hexagonal pattern of air holes, except at the center. must be solved with vector solver. """
    def __init__(self,rhole,rclad,nclad,spacing,hole_res,clad_res,hole_mesh_size=None,clad_mesh_size=None,nhole=1.,rcore=0):
        """
        ARGS:
            rhole (float): the radius of each air hole perforating the fiber
            rclad (float): the outer cladding radius of the fiber
            nclad (float): the index of the cladding material
            spacing (float): the physical spacing between holes
            hole_res (int): the # of line segments used to resolve each hole boundary
            clad_res (int): the # of line segments used to resolve the outer cladding boundary
            hole_mesh_size (float): target mesh size inside holes
            clad_mesh_size (float): target mesh size inside cladding, but outside holes
            nhole (float): index of the holes, default 1.
            rcore (float): the "core radius" of the fiber. holes whose centers are within this radius from the origin are not generated. default 0 (no central hole).
        """    
        
        # get air hole positions
        layers = int(rclad/spacing)
        xa = ya = np.linspace(-layers*spacing,layers*spacing,2*layers+1,endpoint=True)
        xg , yg = np.meshgrid(xa,ya)

        yg *= np.sqrt(3)/2
        if layers%2==1:
            xg[::2, :] += 0.5 * spacing
        else:
            xg[1::2, :] += 0.5 * spacing

        rg = np.sqrt(xg*xg + yg*yg)
        xhole , yhole = xg[rg < rclad-rhole].flatten() , yg[rg < rclad-rhole].flatten()

        # make holes
        holes = []
        for xh,yh in zip(xhole,yhole):
            if xh*xh+yh*yh <= rcore*rcore:
                continue
            hole = Circle(nhole,None)
            hole.make_points(rhole,hole_res,(xh,yh))
            hole.mesh_size = hole_mesh_size
            holes.append(hole)

        # make cladding
        cladding = Circle(nclad,"cladding")
        cladding.make_points(rclad,clad_res)
        cladding.mesh_size = clad_mesh_size

        super().__init__([cladding,Prim2DArray(holes,"holes")])

class PhotonicBandgapFiber(Waveguide):
    """ an optical fiber filled with a hexagonal pattern of air holes, except at the center. must be solved with vector solver. """
    def __init__(self,rvoid,rhole,rclad,nclad,spacing,hole_res,clad_res,hole_mesh_size=None,clad_mesh_size=None,nhole=1.):
        """
        ARGS:
            rvoid (float): the radius of the central air hole
            rhole (float): the radius of the cladding air holes perforating the fiber
            rclad (float): the outer cladding radius of the fiber
            nclad (float): the index of the cladding material
            spacing (float): the physical spacing between holes
            hole_res (int): the # of line segments used to resolve each hole boundary
            clad_res (int): the # of line segments used to resolve the outer cladding boundary
            hole_mesh_size (float): target mesh size inside holes
            clad_mesh_size (float): target mesh size inside cladding, but outside holes
            nhole (float): index of the holes, default 1.
        """    
        
        # get air hole positions
        layers = int(rclad/spacing)+1
        xa = ya = np.linspace(-layers*spacing,layers*spacing,2*layers+1,endpoint=True)
        xg , yg = np.meshgrid(xa,ya)

        yg *= np.sqrt(3)/2
        if layers%2==1:
            xg[::2, :] += 0.5 * spacing
        else:
            xg[1::2, :] += 0.5 * spacing

        rg = np.sqrt(xg*xg + yg*yg)
        xhole , yhole = xg[rg < rclad-rhole].flatten() , yg[rg < rclad-rhole].flatten()

        # make holes
        holes = []
        overlapped_holes = []
        for xh,yh in zip(xhole,yhole):
            hole = Circle(nhole,None)
            hole.make_points(rhole,hole_res,(xh,yh))
            hole.mesh_size = hole_mesh_size
            if np.sqrt(xh*xh+yh*yh)-rhole <= rvoid:
                overlapped_holes.append(hole)
            else:
                holes.append(hole)

        # make cladding
        cladding = Circle(nclad,"cladding")
        cladding.make_points(rclad,clad_res)
        cladding.mesh_size = clad_mesh_size
        
        # make center void
        center = Circle(nhole,None)
        center.make_points(rvoid,int(rvoid/rhole)*hole_res)
        overlapped_holes.append(center)

        void = Prim2DUnion(overlapped_holes,"void")
        void.mesh_size = hole_mesh_size

        hole_array = Prim2DArray(holes,"holes")
        hole_array.mesh_size = hole_mesh_size

        super().__init__([cladding,[hole_array,void]])


class FiberBundleLantern(Waveguide):
    """Photonic lantern with hexagonal arrangement of individual fibers - WaveSolve compatible"""

    def __init__(self, r_jack, r_fiber_clad, r_core, n_rings, n_core, n_clad,
                 core_res=16, clad_res=32, jack_res=None,
                 spacing_factor=2.0, include_center=True,
                 taper_ratio=1.0, r_target_mmcore_size=None,
                 core_mesh_size=None, clad_mesh_size=None,
                 n_jack=None, center_clad_factor=1.5,
                 ring_clad_factors=None):
        """
        Initialize hexagonal fiber bundle photonic lantern.

        Args:
            r_jack: jacket radius
            r_fiber_clad: individual fiber cladding radius
            r_core: core radius (final size after taper)
            n_rings: number of hexagonal rings
            n_core: core refractive index
            n_clad: cladding refractive index
            core_res: resolution for each core circle
            clad_res: resolution for each fiber cladding circle
            jack_res: resolution for jacket circle (default clad_res/2)
            spacing_factor: multiplier for fiber spacing (center-to-center)
            include_center: whether to include center fiber
            taper_ratio: scaling factor (initial_size/final_size)
            r_target_mmcore_size: desired MM core size. Will override taper ratio.
            core_mesh_size: target mesh size in cores
            clad_mesh_size: target mesh size in cladding
            n_jack: jacket refractive index (default same as cladding)
            center_clad_factor: factor to enlarge center fiber cladding (simulates fusing)
            ring_clad_factors: list/dict of cladding scaling factors for each ring
                              Can be:
                              - list: [ring0_factor, ring1_factor, ring2_factor, ...]
                              - dict: {0: ring0_factor, 1: ring1_factor, ...}
                              - None: use center_clad_factor for ring 0, 1.0 for others
        """
        if jack_res is None:
            jack_res = int(clad_res / 2)
        if n_jack is None:
            n_jack = n_clad

        # Process ring cladding factors
        self.ring_clad_factors = self._process_ring_clad_factors(
            ring_clad_factors, n_rings, center_clad_factor, include_center
        )

        # Calculate taper ratio based on target bundle size
        if r_target_mmcore_size is not None:
            # Calculate the radius of the outermost fiber bundle without taper
            original_bundle_radius = self._calculate_bundle_radius(
                n_rings, r_fiber_clad, spacing_factor, include_center
            )
            taper_ratio = r_target_mmcore_size / original_bundle_radius

        # Apply taper ratio to all dimensions
        r_fiber_clad_tapered = r_fiber_clad * taper_ratio
        r_core_tapered = r_core * taper_ratio
        r_jack_tapered = r_jack * taper_ratio

        # Calculate fiber spacing (center-to-center distance)
        spacing = spacing_factor * r_fiber_clad_tapered
        fiber_positions, fiber_rings = self._hex_grid_positions_with_rings(
            n_rings, spacing, include_center
        )

        # Create jacket
        jacket = Circle(n_jack, "jacket")
        jacket.make_points(r_jack_tapered, jack_res)
        jacket.mesh_size = clad_mesh_size

        # Create individual fiber claddings and cores
        fiber_claddings = []
        cores = []

        for i, (pos, ring_idx) in enumerate(zip(fiber_positions, fiber_rings)):
            # Get cladding scaling factor for this ring
            clad_factor = self.ring_clad_factors.get(ring_idx, 1.0)
            clad_radius = r_fiber_clad_tapered * clad_factor

            # Create fiber cladding
            fiber_clad = Circle(n_clad, "cladding")
            fiber_clad.make_points(clad_radius, clad_res, center=pos)
            fiber_clad.mesh_size = clad_mesh_size
            fiber_claddings.append(fiber_clad)

            # Create core at the same position
            core = Circle(n_core, "core")
            core.make_points(r_core_tapered, core_res, center=pos)
            core.mesh_size = core_mesh_size
            cores.append(core)

        # Create arrays for claddings and cores
        fiber_clad_array_Union = Prim2DUnion(fiber_claddings, "cladding")
        fiber_clad_array_Union.mesh_size = clad_mesh_size

        core_array = Prim2DArray(cores, "core")

        # Store metadata
        self.n_fibers = len(fiber_positions)
        self.fiber_positions = fiber_positions
        self.fiber_rings = fiber_rings
        self.taper_ratio = taper_ratio
        self.spacing = spacing
        self.r_fiber_clad = r_fiber_clad_tapered
        self.center_clad_factor = center_clad_factor
        self.bundle_radius = self._calculate_actual_bundle_radius(fiber_positions, r_fiber_clad_tapered)

        # Initialize waveguide with layers (jacket, fiber claddings, cores)
        super().__init__([jacket, fiber_clad_array_Union, core_array])

    def _process_ring_clad_factors(self, ring_clad_factors, n_rings, center_clad_factor, include_center):
        """Process and validate ring cladding factors"""
        factors = {}

        if ring_clad_factors is None:
            # Default behavior: center gets center_clad_factor, others get 1.0
            if include_center:
                factors[0] = center_clad_factor
            for ring in range(1, n_rings + 1):
                factors[ring] = 1.0

        elif isinstance(ring_clad_factors, (list, tuple)):
            # List format: [ring0, ring1, ring2, ...]
            start_ring = 0 if include_center else 1
            for i, factor in enumerate(ring_clad_factors):
                ring_idx = start_ring + i
                if ring_idx <= n_rings:
                    factors[ring_idx] = factor

            # Fill in missing rings with 1.0
            for ring in range(start_ring, n_rings + 1):
                if ring not in factors:
                    factors[ring] = 1.0

        elif isinstance(ring_clad_factors, dict):
            # Dictionary format: {ring_idx: factor}
            factors = ring_clad_factors.copy()

            # Fill in missing rings with 1.0
            start_ring = 0 if include_center else 1
            for ring in range(start_ring, n_rings + 1):
                if ring not in factors:
                    factors[ring] = 1.0

        else:
            raise ValueError("ring_clad_factors must be None, list, tuple, or dict")

        return factors

    def _hex_grid_positions_with_rings(self, n_rings, spacing, include_center=True):
        """Generate hexagonal grid positions for fiber centers with ring information"""
        positions = []
        rings = []

        if include_center:
            positions.append((0, 0))
            rings.append(0)

        for ring in range(1, n_rings + 1):
            for i in range(6 * ring):
                angle = 2 * np.pi * i / (6 * ring)
                edge = int(i / ring)
                pos_on_edge = i % ring
                edge_angle = edge * np.pi / 3
                edge_dir = (edge + 2) * np.pi / 3

                x = ring * spacing * np.cos(edge_angle) + pos_on_edge * spacing * np.cos(edge_dir)
                y = ring * spacing * np.sin(edge_angle) + pos_on_edge * spacing * np.sin(edge_dir)

                positions.append((x, y))
                rings.append(ring)

        return positions, rings

    def _calculate_bundle_radius(self, n_rings, r_fiber_clad, spacing_factor, include_center=True):
        """Calculate the radius of the enscribing circle for the fiber bundle"""
        if n_rings == 0:
            return r_fiber_clad

        # Calculate fiber spacing
        spacing = spacing_factor * r_fiber_clad

        # Distance from center to outermost fiber centers
        outermost_distance = n_rings * spacing

        # Add the fiber cladding radius to get the enscribing circle
        bundle_radius = outermost_distance + r_fiber_clad

        return bundle_radius

    def _calculate_actual_bundle_radius(self, fiber_positions, r_fiber_clad_tapered):
        """Calculate the actual bundle radius from fiber positions"""
        if not fiber_positions:
            return r_fiber_clad_tapered

        # Find the maximum distance from center to any fiber edge
        max_distance = 0
        for pos in fiber_positions:
            fiber_center_distance = np.sqrt(pos[0] ** 2 + pos[1] ** 2)
            fiber_edge_distance = fiber_center_distance + r_fiber_clad_tapered
            max_distance = max(max_distance, fiber_edge_distance)

        return max_distance

    def _hex_grid_positions(self, n_rings, spacing, include_center=True):
        """Generate hexagonal grid positions for fiber centers (legacy method)"""
        positions, _ = self._hex_grid_positions_with_rings(n_rings, spacing, include_center)
        return positions

    def make_mesh(self, algo=6, order=2, adaptive=True):
        """Generate mesh with enhanced control for fiber bundle lanterns"""
        mesh = super().make_mesh(algo, order, adaptive)

        # Add fiber bundle-specific metadata
        mesh.field_data.update({
            "n_fibers": self.n_fibers,
            "taper_ratio": self.taper_ratio,
            "spacing": self.spacing,
            "center_clad_factor": self.center_clad_factor,
            "ring_clad_factors": self.ring_clad_factors,
            "bundle_radius": self.bundle_radius
        })

        return mesh

    def get_fiber_info(self):
        """Return information about individual fibers"""
        info = {
            'n_fibers': self.n_fibers,
            'fiber_positions': self.fiber_positions,
            'fiber_rings': self.fiber_rings,
            'fiber_cladding_radius': self.r_fiber_clad,
            'center_enlarged': self.center_clad_factor > 1.0,
            'center_clad_factor': self.center_clad_factor,
            'ring_clad_factors': self.ring_clad_factors,
            'bundle_radius': self.bundle_radius
        }
        return info

    def get_ring_info(self):
        """Return information about rings and their cladding factors"""
        ring_info = {}
        for ring_idx, factor in self.ring_clad_factors.items():
            fiber_count = 1 if ring_idx == 0 else 6 * ring_idx
            ring_info[ring_idx] = {
                'cladding_factor': factor,
                'fiber_count': fiber_count,
                'effective_clad_radius': self.r_fiber_clad * factor
            }
        return ring_info

class MCFPhotonicLantern(Waveguide):
    """Photonic lantern with hexagonal grid of core, i.e. MCF based lantern."""

    def __init__(self, r_jack, r_clad, r_core, n_rings, n_core, n_clad,
                 core_res=16, clad_res=32, jack_res=None,
                 spacing_factor=2.2, include_center=True,
                 taper_ratio=1.0, r_target_cladding_size=None,
                 core_mesh_size=None, clad_mesh_size=None,
                 n_jack=None):
        """
        Initialize MCF photonic lantern.

        Args:
            r_jack: jacket radius
            r_clad: cladding radius of MCF
            r_core: core radius (final size after taper)
            n_rings: number of hexagonal rings
            n_core: core refractive index
            n_clad: cladding refractive index
            core_res: resolution for each core circle
            clad_res: resolution for cladding circle
            jack_res: resolution for jacket circle (default clad_res/2)
            spacing_factor: multiplier for core spacing
            include_center: whether to include center core
            taper_ratio: scaling factor (initial_size/final_size)
            r_target_cladding_size: desired MM core size. Will overide taper ratio.
            core_mesh_size: target mesh size in cores
            clad_mesh_size: target mesh size in cladding
            n_jack: jacket refractive index (default same as cladding)
        """
        if jack_res is None:
            jack_res = int(clad_res/2)
        if n_jack is None:
            n_jack = n_clad

        # Calculate core positions based on taper
        if r_target_cladding_size is not None:
            taper_ratio = r_target_cladding_size/r_clad

        spacing_base = r_core if taper_ratio == 1.0 else r_core * taper_ratio
        r_clad=r_clad*taper_ratio
        r_core=r_core*taper_ratio
        r_jack=r_jack*taper_ratio

        spacing = spacing_factor * spacing_base
        core_positions = self._hex_grid_positions(n_rings, spacing, include_center)

        # Create jacket
        jacket = Circle(n_jack, "jacket")
        jacket.make_points(r_jack, jack_res)
        jacket.mesh_size = clad_mesh_size

        # Create cladding
        cladding = Circle(n_clad, "cladding")
        cladding.make_points(r_clad, clad_res)
        cladding.mesh_size = clad_mesh_size

        # Create cores
        cores = []
        for i, pos in enumerate(core_positions):
            core = Circle(n_core, None)  # Individual label added later
            core.make_points(r_core, core_res, center=pos)
            core.mesh_size = core_mesh_size
            cores.append(core)

        # Create core array
        core_array = Prim2DArray(cores, "cores")

        # Store metadata
        self.n_cores = len(cores)
        self.core_positions = core_positions
        self.taper_ratio = taper_ratio
        self.spacing = spacing

        # Initialize waveguide with layers
        super().__init__([jacket, cladding, core_array])

        self.min_mesh_size = min(r_core/4, 0.1)
        self.max_mesh_size = r_clad/5

    def _hex_grid_positions(self, n_rings, spacing, include_center=True):
        """Generate hexagonal grid positions"""
        positions = []

        if include_center:
            positions.append((0, 0))

        for ring in range(1, n_rings + 1):
            for i in range(6 * ring):
                angle = 2 * np.pi * i / (6 * ring)
                edge = int(i / ring)
                pos_on_edge = i % ring
                edge_angle = edge * np.pi / 3
                edge_dir = (edge + 2) * np.pi / 3

                x = ring * spacing * np.cos(edge_angle) + pos_on_edge * spacing * np.cos(edge_dir)
                y = ring * spacing * np.sin(edge_angle) + pos_on_edge * spacing * np.sin(edge_dir)

                positions.append((x, y))

        return positions

    def make_mesh(self, algo=6, order=2, adaptive=True):
        """Generate mesh with enhanced control for photonic lanterns"""
        mesh = super().make_mesh(algo, order, adaptive)

        # Add lantern-specific metadata
        mesh.field_data.update({
            "n_cores": self.n_cores,
            "taper_ratio": self.taper_ratio,
            "spacing": self.spacing
        })

        return mesh
#endregion