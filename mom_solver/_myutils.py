# -*- coding: utf-8 -*-
import numpy as np
import numpy.linalg
import scipy
import scipy.sparse
import scipy.sparse.linalg
import re
from scipy.sparse.linalg import LinearOperator
from scipy.sparse import coo_matrix

'''
程序准备部分
'''
# 三角形的高斯积分
class Cubature(object):
    def __init__(self, degree, rule): # 初始化
    	self.point_ = []
    	self.weight_ = []
    	self.degree_ = degree
    	self.rule_ = rule
    	self.numPoint_ = 0
    	self.currentPointIndex_ = 0
    	self.setupTriangleCubtature()
    	pass
    def point(self, vertex,pointIndex): # 输出三角形的高斯点
        return self.point_[pointIndex][0] * np.array(vertex[0]) + self.point_[pointIndex][1] *  np.array(vertex[1]) + self.point_[pointIndex][2] *  np.array(vertex[2])   
        pass
    def weight(self, area, pointIndex):  # 输出三角形高斯积分的权重
        return area*2. * self.weight_[pointIndex];
        pass
    def kesi(self, pointIndex, i): # 标准三角形的面积坐标下高斯点位置
        temp = {0:self.point_[pointIndex][0],1:self.point_[pointIndex][1],2:self.point_[pointIndex][2]}
        return temp[i]
        pass
    def numPoint(self): # 输出高斯点的个数
        return self.numPoint_
        pass
    def degree(self): # 输出高斯积分的阶次
        return self.degree_
        pass
    def rule(self): # 输出高斯积分的方法
        return self.rule_
        pass
    def setupComplete(self):
        pass
    def setupTriangleCubtature(self): # 初始化用的内部函数
        self.findNumPointInTriangleCubature()
        self.point_ = [(0,0,0) for x in xrange(self.numPoint_)]
        self.weight_ = [0 for x in xrange(self.numPoint_)]
        self.findPointAndWeightInTriangle()
        pass
    def findNumPointInTriangleCubature(self): # 计算高斯点的个数
        numPoints = {101:1,201:3,301:4,302:6,401:6,501:7,601:12,602:12,701:12,801:16,901:19,1001:25,1002:25,1101:28}
        try:
            self.numPoint_ = numPoints[100*self.degree_ + self.rule_]
        except:
            print "Neither degree or rule is correct!"
            raise()
        pass
                
    def findPointAndWeightInTriangle(self): # 标准三角形的高斯点和权重
        if (100*self.degree_ + self.rule_) == 101:
            self.SIMPLEX2_FULLY1(0.33333333333333333,0.33333333333333333, 0.50000000000000000)
            self.currentPointIndex_ = 0
        elif (100*self.degree_ + self.rule_) == 201:
            self.SIMPLEX2_FULLY3(0.16666666666666666,0.16666666666666666, 0.16666666666666666)
            self.currentPointIndex_ = 0
        elif (100*self.degree_ + self.rule_) ==  301:
            self.SIMPLEX2_FULLY1(0.33333333333333333,0.33333333333333333, -0.28125000000000000)
            self.SIMPLEX2_FULLY3(0.20000000000000000,0.20000000000000000,  0.26041666666666667)
            self.currentPointIndex_ = 0
        elif (100*self.degree_ + self.rule_) ==  302:
            self.SIMPLEX2_FULLY6(0.10903900907287721,0.23193336855303057, 0.083333333333333333)
        elif (100*self.degree_ + self.rule_) ==  401:
            self.SIMPLEX2_FULLY3(0.091576213509770743,0.091576213509770743, 0.054975871827660933)
            self.SIMPLEX2_FULLY3(0.44594849091596488, 0.44594849091596488,  0.11169079483900573)
            self.currentPointIndex_ = 0
        elif (100*self.degree_ + self.rule_) ==  501:
            self.SIMPLEX2_FULLY1(0.33333333333333333,0.33333333333333333, 0.11250000000000000)
            self.SIMPLEX2_FULLY3(0.10128650732345633,0.10128650732345633, 0.062969590272413576)
            self.SIMPLEX2_FULLY3(0.47014206410511508,0.47014206410511508, 0.066197076394253090)
            self.currentPointIndex_ = 0
        elif (100*self.degree_ + self.rule_) ==  601:
            self.SIMPLEX2_FULLY3(0.063089014491502228,0.063089014491502228, 0.025422453185103408)
            self.SIMPLEX2_FULLY3(0.24928674517091042, 0.24928674517091042,  0.058393137863189683)
            self.SIMPLEX2_FULLY6(0.053145049844816947,0.31035245103378440,  0.041425537809186787)
            self.currentPointIndex_ = 0
        elif (100*self.degree_ + self.rule_) ==  602:
            self.SIMPLEX2_FULLY3(0.21942998254978296,0.21942998254978296, 0.085666562076490515)
            self.SIMPLEX2_FULLY3(0.48013796411221504,0.48013796411221504, 0.040365544796515489)
            self.SIMPLEX2_FULLY6(0.83900925971479105,0.14161901592396815, 0.020317279896830331)
            self.currentPointIndex_ = 0
        elif (100*self.degree_ + self.rule_) ==  701:
            self.SIMPLEX2_ROTATIONAL3(0.062382265094402118,0.067517867073916085, 0.026517028157436251)
            self.SIMPLEX2_ROTATIONAL3(0.055225456656926611,0.32150249385198182,  0.043881408714446055)
            self.SIMPLEX2_ROTATIONAL3(0.034324302945097146,0.66094919618673565,  0.028775042784981585)
            self.SIMPLEX2_ROTATIONAL3(0.51584233435359177,0.27771616697639178, 0.067493187009802774)
            self.currentPointIndex_ = 0
        elif (100*self.degree_ + self.rule_) ==  801:
            self.SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.072157803838893584)
            self.SIMPLEX2_FULLY3(0.17056930775176020, 0.17056930775176020,  0.051608685267359125)
            self.SIMPLEX2_FULLY3(0.050547228317030975,0.050547228317030975, 0.016229248811599040)
            self.SIMPLEX2_FULLY3(0.45929258829272315, 0.45929258829272315,  0.047545817133642312)
            self.SIMPLEX2_FULLY6(0.72849239295540428, 0.26311282963463811,  0.013615157087217497)
            self.currentPointIndex_ = 0
        elif (100*self.degree_ + self.rule_) ==  901:
            self.SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.048567898141399416)
            self.SIMPLEX2_FULLY3(0.48968251919873762, 0.48968251919873762,  0.015667350113569535)
            self.SIMPLEX2_FULLY3(0.43708959149293663, 0.43708959149293663,  0.038913770502387139)
            self.SIMPLEX2_FULLY3(0.18820353561903273, 0.18820353561903273,  0.039823869463605126)
            self.SIMPLEX2_FULLY3(0.044729513394452709,0.044729513394452709, 0.012788837829349015)
            self.SIMPLEX2_FULLY6(0.74119859878449802, 0.036838412054736283, 0.021641769688644688)
            self.currentPointIndex_ = 0
        elif (100*self.degree_ + self.rule_) ==  1001:
            self.SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.039947252370619853)
            self.SIMPLEX2_FULLY3(0.42508621060209057, 0.42508621060209057,  0.035561901116188667)
            self.SIMPLEX2_FULLY3(0.023308867510000190,0.023308867510000190, 4.1119093452320977e-3)
            self.SIMPLEX2_FULLY6(0.62830740021349255, 0.22376697357697300,  0.022715296148085009)
            self.SIMPLEX2_FULLY6(0.61131382618139764, 0.35874014186443146,  0.018679928117152638)
            self.SIMPLEX2_FULLY6(0.82107206998562937, 0.14329537042686714,  0.015443328442281994)
            self.currentPointIndex_ = 0
        elif (100*self.degree_ + self.rule_) ==  1002:
            self.SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.040871664573142983)
            self.SIMPLEX2_FULLY3(0.14216110105656438, 0.14216110105656438,  0.022978981802372364)
            self.SIMPLEX2_FULLY3(0.032055373216943512,0.032055373216943512, 6.6764844065747831e-3)
            self.SIMPLEX2_FULLY6(0.53005411892734402, 0.32181299528883542,  0.031952453198212022)
            self.SIMPLEX2_FULLY6(0.60123332868345924, 0.36914678182781098,  0.017092324081479714)
            self.SIMPLEX2_FULLY6(0.80793060092287906, 0.16370173373718249,  0.012648878853644192)
            self.currentPointIndex_ = 0
        elif (100*self.degree_ + self.rule_) ==  1101:
            self.SIMPLEX2_FULLY6(0.85887028128263670, 0.14112971871736329,  3.6811918916502771e-3)
            self.SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.043988650581116119)
            self.SIMPLEX2_FULLY3(0.025989140928287395,0.025989140928287395, 4.3721557768680115e-3)
            self.SIMPLEX2_FULLY3(0.094287502647922495,0.094287502647922495, 0.019040785996967468)
            self.SIMPLEX2_FULLY3(0.49463677501721381, 0.49463677501721381,  9.4277240280656460e-3) 
            self.SIMPLEX2_FULLY3(0.20734338261451133, 0.20734338261451133,  0.036079848772369763)
            self.SIMPLEX2_FULLY3(0.43890780570049209, 0.43890780570049209,  0.034664569352767949)
            self.SIMPLEX2_FULLY6(0.67793765488259040, 0.044841677589130443, 0.020528157714644283)
            self.currentPointIndex_ = 0
        else:
            self.SIMPLEX2_FULLY6(0.85887028128263670, 0.14112971871736329,  3.6811918916502771e-3)
            self.SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.043988650581116119)
            self.SIMPLEX2_FULLY3(0.025989140928287395,0.025989140928287395, 4.3721557768680115e-3)
            self.SIMPLEX2_FULLY3(0.094287502647922495,0.094287502647922495, 0.019040785996967468)
            self.SIMPLEX2_FULLY3(0.49463677501721381, 0.49463677501721381,  9.4277240280656460e-3) 
            self.SIMPLEX2_FULLY3(0.20734338261451133, 0.20734338261451133,  0.036079848772369763)
            self.SIMPLEX2_FULLY3(0.43890780570049209, 0.43890780570049209,  0.034664569352767949)
            self.SIMPLEX2_FULLY6(0.67793765488259040, 0.044841677589130443, 0.020528157714644283)
            self.currentPointIndex_ = 0
        pass
    def SIMPLEX2_FULLY1(self, x,y, weight): # 内部函数
        coord1 = x
        coord2 = y
        coord3 = 1.0 - coord1 - coord2
        self.point_[self.currentPointIndex_] = (coord1, coord2, coord3)
        self.weight_[self.currentPointIndex_] = weight
        self.currentPointIndex_ = self.currentPointIndex_+1
        pass
    def SIMPLEX2_FULLY3(self, x,y, weight): # 内部函数
        coord1 = x
        coord2 = y
        coord3 = 1.0 - coord1 - coord2
        self.point_[self.currentPointIndex_] = (coord1, coord2, coord3)
        self.weight_[self.currentPointIndex_] = weight
        self.currentPointIndex_ = self.currentPointIndex_ + 1
        self.point_[self.currentPointIndex_] = (coord3, coord1, coord2)
        self.weight_[self.currentPointIndex_] = weight
        self.currentPointIndex_ = self.currentPointIndex_ + 1
        self.point_[self.currentPointIndex_] = (coord2, coord3, coord1)
        self.weight_[self.currentPointIndex_] = weight
        self.currentPointIndex_ = self.currentPointIndex_ + 1
        pass
    def SIMPLEX2_FULLY6(self, x,y, weight): # 内部函数
        coord1 = x
        coord2 = y
        coord3 = 1.0 - coord1 - coord2
        self.point_[self.currentPointIndex_] = (coord1, coord2, coord3)
        self.weight_[self.currentPointIndex_] = weight
        self.currentPointIndex_ = self.currentPointIndex_ + 1

        self.point_[self.currentPointIndex_] = (coord1, coord3, coord2) 
        self.weight_[self.currentPointIndex_] = weight
        self.currentPointIndex_ = self.currentPointIndex_ + 1

        self.point_[self.currentPointIndex_] = (coord2, coord1, coord3) 
        self.weight_[self.currentPointIndex_] = weight
        self.currentPointIndex_ = self.currentPointIndex_ + 1

        self.point_[self.currentPointIndex_] = (coord3, coord1, coord2) 
        self.weight_[self.currentPointIndex_] = weight
        self.currentPointIndex_ = self.currentPointIndex_ + 1

        self.point_[self.currentPointIndex_] = (coord2, coord3, coord1) 
        self.weight_[self.currentPointIndex_] = weight
        self.currentPointIndex_ = self.currentPointIndex_ + 1

        self.point_[self.currentPointIndex_] = (coord3, coord2, coord1) 
        self.weight_[self.currentPointIndex_] = weight
        self.currentPointIndex_ = self.currentPointIndex_ + 1
        pass
    def SIMPLEX2_ROTATIONAL3(self, x,y, weight): # 内部函数
        coord1 = x
        coord2 = y
        coord3 = 1.0 - coord1 - coord2;
        self.point_[self.currentPointIndex_] = (coord1, coord2, coord3); 
        self.weight_[self.currentPointIndex_] = weight;
        self.currentPointIndex_ = self.currentPointIndex_ + 1
        self.point_[self.currentPointIndex_] = (coord3, coord1, coord2); 
        self.weight_[self.currentPointIndex_] = weight;
        self.currentPointIndex_ = self.currentPointIndex_ + 1
        self.point_[self.currentPointIndex_] = (coord2, coord3, coord1); 
        self.weight_[self.currentPointIndex_] = weight;
        self.currentPointIndex_ = self.currentPointIndex_ + 1
        pass

# 定义三角形方法
class Triangle(object):
    def __init__(self,triPoints):
        self.triPoints = triPoints
    def test(self,r): # 判断是否在三角型内部
        pointA = self.triPoints[0]
        pointB = self.triPoints[1]
        pointC = self.triPoints[2]
        area1 = self.area(pointA,pointB,r)
        area2 = self.area(pointB,pointC,r)
        area3 = self.area(pointC,pointA,r)
        areaT = self.area(pointA,pointB,pointC)
        if np.sum([area1,area2,area3]) <= areaT*1.0001:
            return True
        elif np.sum([area1,area2,area3]) > areaT*1.0001:
            return False
        else:
            raise
    def test2(self,r):
        pointA = np.array(self.triPoints[0])
        pointB = np.array(self.triPoints[1])
        pointC = np.array(self.triPoints[2])
        center = (pointA + pointB + pointC)/3.
        threshold = max([np.linalg.norm(np.array(self.triPoints[ii])-center,2) for ii in xrange(3)])
        if np.linalg.norm(center-np.array(r),2)<=3.*threshold:
            return True
        else:
            return False

    def area(self,A,B,C): # 计算三角型面积
        v1 = np.array(A)-np.array(B)
        v2 = np.array(C)-np.array(B)
        return np.linalg.norm(np.cross(v1,v2))/2.0

# 奇异性处理，参考wilton TAP 1983的论文
class SiguralityInt(object):
    def __init__(self):
        pass


    def _R_Polygon(self, polysurf, r):
        try:
            norm_ = np.cross(polysurf[1]-polysurf[0], polysurf[2]-polysurf[0])
            norm_ = norm_/np.linalg.norm(norm_,2)
            l_ = []
            u_ = []
            r_ = []
            R = []
            for ii in xrange(3):
                temp = polysurf[(ii+1)%3]-polysurf[ii]
                l_.append( temp/np.linalg.norm(temp,2) )
                u_.append( np.cross(l_[ii],norm_) )
                r_.append( polysurf[ii]-r )
                R.append( np.linalg.norm(r_[ii],2) )
            d = r_[0].dot(norm_)
            P_ = []
            for ii in xrange(3):
                temp = d*norm_
                P_.append (r_[ii]-temp)
            P0 = []
            lpos = []
            lneg = []
            R0 = []
            for ii in xrange(3):
                P0.append(abs(np.dot(P_[ii], u_[ii])))
                lpos.append (np.dot(P_[(ii+1)%3], l_[ii]))
                lneg.append( np.dot(P_[ii], l_[ii]) )
                R0.append( scipy.sqrt(P0[ii]*P0[ii]+d*d) )
            result =0
            for ii in xrange(3):
                noise = 1.e-10*np.linalg.norm(l_[ii],2)
                if (R0[ii]<noise):  
                    continue
                if (np.dot(P_[ii], u_[ii])>0):
                    sign = 1
                else:
                    sign = -1
                lg=0.0
                if (R[ii]+lneg[ii]>noise):
                    lg=	scipy.log(R[(ii+1)%3]+lpos[ii]) - scipy.log(R[ii]+lneg[ii])
                Im1=lg
                Ip1=0.5*(R[(ii+1)%3]*lpos[ii]-R[ii]*lneg[ii]+R0[ii]*R0[ii]*Im1)
                if abs(d)>noise:
                    result = result + sign*(P0[ii]*Ip1+abs(d)*abs(d)*(P0[ii]*Im1-abs(d)*(scipy.arctan(P0[ii]*lpos[ii]/(R0[ii]*R0[ii]+abs(d)*R[(ii+1)%3]))-scipy.arctan(P0[ii]*lneg[ii]/(R0[ii]*R0[ii]+abs(d)*R[ii])))))/3.0
                else:
                    result = result + sign*(P0[ii]*Ip1)/3.0
            return result
            pass
        except:
            raise

    # 奇异核的解析积分
    def _1_R_Polygon(self, polysurf, r): # 参见经典文献中的公式
        try:
            norm_ = np.cross(polysurf[1]-polysurf[0], polysurf[2]-polysurf[0])
            norm_ = norm_/np.linalg.norm(norm_,2)
            l_ = []
            u_ = []
            r_ = []
            R = []
            for ii in xrange(3):
                temp = polysurf[(ii+1)%3]-polysurf[ii]
                l_.append( temp/np.linalg.norm(temp,2) )
                u_.append( np.cross(l_[ii],norm_) )
                r_.append( polysurf[ii]-r )
                R.append( np.linalg.norm(r_[ii],2) )
            d = r_[0].dot(norm_)
            P_ = []
            for ii in xrange(3):
                temp = d*norm_
                P_.append (r_[ii]-temp)
            P0 = []
            lpos = []
            lneg = []
            R0 = []
            for ii in xrange(3):
                P0.append(abs(np.dot(P_[ii], u_[ii])))
                lpos.append (np.dot(P_[(ii+1)%3], l_[ii]))
                lneg.append( np.dot(P_[ii], l_[ii]) )
                R0.append( scipy.sqrt(P0[ii]*P0[ii]+d*d) )

            result =0
            for ii in xrange(3):
                noise = 1.e-10*np.linalg.norm(l_[ii],2)
                if (R0[ii]<noise):  
                    continue
                if (np.dot(P_[ii], u_[ii])>0):
                    sign = 1
                else:
                    sign = -1
                lg=0.0
                if (R[ii]+lneg[ii]>noise):
                    lg=	scipy.log(R[(ii+1)%3]+lpos[ii]) - scipy.log(R[ii]+lneg[ii])
                if abs(d)>noise:
                    result = result + sign*( P0[ii]*lg - abs(d)*( scipy.arctan(P0[ii]*lpos[ii]/(R0[ii]*R0[ii]+abs(d)*R[(ii+1)%3]))-scipy.arctan(P0[ii]*lneg[ii]/(R0[ii]*R0[ii]+abs(d)*R[ii])) ) )
                else:
                    result = result + sign*P0[ii]*lg
            return result
        except:
            raise
    def _grad_R_Polygon(self, polysurf, r): # 参见经典文献中的公式
        try:
            norm_ = np.cross(polysurf[1]-polysurf[0], polysurf[2]-polysurf[0])
            norm_ = norm_/np.linalg.norm(norm_,2)
            l_ = []
            u_ = []
            r_ = []
            R = []
            for ii in xrange(3):
                temp = polysurf[(ii+1)%3]-polysurf[ii]
                l_.append( temp/np.linalg.norm(temp,2) )
                u_.append( np.cross(l_[ii],norm_) )
                r_.append( polysurf[ii]-r )
                R.append( np.linalg.norm(r_[ii],2) )
            d = r_[0].dot(norm_)
            P_ = []

            for ii in xrange(3):
                temp = d*norm_
                P_.append (r_[ii]-temp)

            P0 = []
            lpos = []
            lneg = []
            R0 = []
            for ii in xrange(3):
                P0.append(abs(np.dot(P_[ii], u_[ii])))
                lpos.append (np.dot(P_[(ii+1)%3], l_[ii]))
                lneg.append( np.dot(P_[ii], l_[ii]) )
                R0.append( scipy.sqrt(P0[ii]*P0[ii]+d*d) )

            result = np.array([0.,0.,0.])
            for jj in xrange(3):
                result[jj] = 0.0
                for ii in xrange(3):
                    noise = 1.e-10*np.linalg.norm(l_[ii],2)
                    if (R0[ii]<noise):  
                        result[jj] = result[jj] + u_[ii][jj]*( R[(ii+1)%3]*lpos[ii]-R[ii]*lneg[ii] )
#                        print "---"
                    else:
                        lg=0.0
                        if (R[ii]+lneg[ii]>noise):
                            lg=	scipy.log(R[(ii+1)%3]+lpos[ii]) - scipy.log(R[ii]+lneg[ii])
#                        print lg
                        result[jj] = result[jj] + u_[ii][jj]*( R0[ii]*R0[ii]*lg+R[(ii+1)%3]*lpos[ii]-R[ii]*lneg[ii] )
#                        print u_[ii][jj]*( R0[ii]*R0[ii]*lg+R[(ii+1)%3]*lpos[ii]-R[ii]*lneg[ii] )
                result[jj] = result[jj]*0.5
                    
            return result
        except:
            raise    

class SiguralityInt_vec(object):
    def __init__(self):
        pass
    def _grad_R_Polygon(self, polysurf, r): # 参见经典文献中的公式
        try:
            norm_ = np.cross(polysurf[:,1,:]-polysurf[:,0,:], polysurf[:,2,:]-polysurf[:,0,:]) # ne*3
            norm_ = norm_/np.sqrt(np.sum(norm_**2, axis=-1).reshape([-1,1]))
            l_ = np.zeros_like(polysurf) # ne*3*3
            u_ = np.zeros_like(polysurf) # ne*3*3
            r_ = np.zeros([r.shape[0],polysurf.shape[0],polysurf.shape[1],polysurf.shape[2]]) # nr*ne*3*3
            R =  np.zeros([r.shape[0],polysurf.shape[0],polysurf.shape[1]]) # nr*ne*3
            for ii in xrange(3):
                temp = polysurf[:,(ii+1)%3,:]-polysurf[:,ii,:]
                l_[:,ii,:] = temp/np.sqrt(np.sum(temp*temp, axis=-1).reshape([-1,1]))
                u_[:,ii,:] =  np.cross(l_[:,ii,:],norm_)
                r_[:,:,ii,:] = polysurf[:,ii,:]-r.reshape([-1,1,3])
                R[:,:,ii] = np.sqrt(np.sum(r_[:,:,ii,:]*r_[:,:,ii,:], axis=-1))
                
            d = np.sum(r_[:,:,0,:]*norm_,axis=-1) # nr*ne
            P_ = np.zeros_like(r_) # nr*ne*3*3
            temp = d.reshape([d.shape[0],d.shape[1],1])*norm_.reshape([1,norm_.shape[0], norm_.shape[1]]) # nr*ne*3
            for ii in xrange(3):
                P_[:,:,ii,:] = r_[:,:,ii,:]-temp

            P0 = np.zeros_like(R) # nr*ne*3
            lpos = np.zeros_like(R)
            lneg = np.zeros_like(R)
            R0 = np.zeros_like(R)
            for ii in xrange(3):
                P0[:,:,ii] = np.abs(np.sum(P_[:,:,ii,:]*u_[:,ii,:],axis=-1))
                lpos[:,:,ii] = np.sum(P_[:,:,(ii+1)%3,:]*l_[:,ii,:],axis=-1)
                lneg[:,:,ii] = np.sum(P_[:,:,ii,:]*l_[:,ii,:],axis=-1)
                R0[:,:,ii] = np.sqrt(P0[:,:,ii]*P0[:,:,ii]+d*d)

            result = np.zeros([r.shape[0],polysurf.shape[0],3])
            for ii in xrange(3):
                noise = 1.e-10*np.sqrt(np.sum(l_[:,ii,:]*l_[:,ii,:],axis=-1))                
                
                temp = R[:,:,(ii+1)%3]*lpos[:,:,ii]-R[:,:,ii]*lneg[:,:,ii]                
                result_branch1 =  u_[:,ii,:].reshape([1,polysurf.shape[0],3])\
                                *temp.reshape([r.shape[0],-1,1])
                
                check2 = (R[:,:,ii]+lneg[:,:,ii])>noise
                lg = np.where(check2,\
                              np.log(R[:,:,(ii+1)%3]+lpos[:,:,ii]) - np.log(R[:,:,ii]+lneg[:,:,ii]),\
                              np.zeros([r.shape[0],polysurf.shape[0]])\
                              )
#                print lg
                temp1111 = R0[:,:,ii]**2*lg\
                            +R[:,:,(ii+1)%3]*lpos[:,:,ii]\
                            -R[:,:,ii]*lneg[:,:,ii]
#                print temp1111
                result_branch2 = u_[:,ii,:].reshape([1,polysurf.shape[0],3])\
                            * temp1111.reshape([r.shape[0],-1,1])
#                print result_branch2
#                
                check1 = (R0[:,:,ii] < noise)
                result = result + np.where(check1.reshape([r.shape[0],polysurf.shape[0],1]), \
                                           result_branch1, \
                                           result_branch2)                                              
            result = result*0.5
            return result
        except Exception as e:
            print e            
            raise
    def _1_R_Polygon(self, polysurf, r): # 参见经典文献中的公式
        try:
            norm_ = np.cross(polysurf[:,1,:]-polysurf[:,0,:], polysurf[:,2,:]-polysurf[:,0,:]) # ne*3
            
            norm_ = norm_/np.sqrt(np.sum(norm_**2, axis=-1).reshape([-1,1]))
            l_ = np.zeros_like(polysurf) # ne*3*3
            u_ = np.zeros_like(polysurf) # ne*3*3
            r_ = np.zeros([r.shape[0],polysurf.shape[0],polysurf.shape[1],polysurf.shape[2]]) # nr*ne*3*3
            R =  np.zeros([r.shape[0],polysurf.shape[0],polysurf.shape[1]]) # nr*ne*3
            for ii in xrange(3):
                temp = polysurf[:,(ii+1)%3,:]-polysurf[:,ii,:]
                l_[:,ii,:] = temp/np.sqrt(np.sum(temp**2, axis=-1).reshape([-1,1]))
                u_[:,ii,:] =  np.cross(l_[:,ii,:],norm_)
                r_[:,:,ii,:] = polysurf[:,ii,:]-r.reshape([-1,1,3])
                R[:,:,ii] = np.sqrt(np.sum(r_[:,:,ii,:]*r_[:,:,ii,:], axis=-1))
           
            d = np.sum(r_[:,:,0,:]*norm_,axis=-1) # nr*ne
            P_ = np.zeros_like(r_) # nr*ne*3*3
            temp = d.reshape([d.shape[0],d.shape[1],1])*norm_.reshape([1,norm_.shape[0], norm_.shape[1]]) # nr*ne*3
            for ii in xrange(3):
                P_[:,:,ii,:] = r_[:,:,ii,:]-temp
            P0 = np.zeros_like(R) # nr*ne*3
            lpos = np.zeros_like(R)
            lneg = np.zeros_like(R)
            R0 = np.zeros_like(R)
            for ii in xrange(3):
                P0[:,:,ii] = np.abs(np.sum(P_[:,:,ii,:]*u_[:,ii,:],axis=-1))
                lpos[:,:,ii] = np.sum(P_[:,:,(ii+1)%3,:]*l_[:,ii,:],axis=-1)
                lneg[:,:,ii] = np.sum(P_[:,:,ii,:]*l_[:,ii,:],axis=-1)
                R0[:,:,ii] = np.sqrt(P0[:,:,ii]*P0[:,:,ii]+d*d)
        
            result = np.zeros([r.shape[0],polysurf.shape[0]])    
            R0__2 = R0**2
            absd = np.abs(d)
            for ii in xrange(3):
                noise = 1.e-10*np.sqrt(np.sum(l_[:,ii,:]**2,axis=-1))                                                
                check2 = (R[:,:,ii]+lneg[:,:,ii])>noise
                lg = np.where(check2,\
                              scipy.log(R[:,:,(ii+1)%3]+lpos[:,:,ii]) \
                              - scipy.log(R[:,:,ii]+lneg[:,:,ii]),\
                              np.zeros([r.shape[0],polysurf.shape[0]])\
                              ) 
                check3 = absd>noise 
                result_branch2 = np.where(check3,\
                                          P0[:,:,ii]*lg - absd*( \
                                            scipy.arctan(P0[:,:,ii]*lpos[:,:,ii]/(R0__2[:,:,ii]+absd*R[:,:,(ii+1)%3]))\
                                            -scipy.arctan(P0[:,:,ii]*lneg[:,:,ii]/(R0__2[:,:,ii]+absd*R[:,:,ii]))), \
                                          P0[:,:,ii]*lg)   
                
                check1 = (R0[:,:,ii] < noise)
                temp_result_add = np.where(check1.reshape([r.shape[0],polysurf.shape[0]]), \
                                           np.zeros_like(result), \
                                           result_branch2)
                sing = np.sum(P_[:,:,ii,:]*u_[:,ii,:], axis=-1)>0
                result = np.where( sing,\
                                  result + temp_result_add,\
                                  result - temp_result_add)

            return result
        except Exception as e:
            print e       
            raise

    def _R_Polygon(self, polysurf, r):
        try:
            norm_ = np.cross(polysurf[:,1,:]-polysurf[:,0,:], polysurf[:,2,:]-polysurf[:,0,:]) # ne*3
            norm_ = norm_/np.sqrt(np.sum(norm_*norm_, axis=-1).reshape([-1,1]))
            l_ = np.zeros_like(polysurf) # ne*3*3
            u_ = np.zeros_like(polysurf) # ne*3*3
            r_ = np.zeros([r.shape[0],polysurf.shape[0],polysurf.shape[1],polysurf.shape[2]]) # nr*ne*3*3
            R =  np.zeros([r.shape[0],polysurf.shape[0],polysurf.shape[1]]) # nr*ne*3
            for ii in xrange(3):
                temp = polysurf[:,(ii+1)%3,:]-polysurf[:,ii,:]
                l_[:,ii,:] = temp/np.sqrt(np.sum(temp*temp, axis=-1).reshape([-1,1]))
                u_[:,ii,:] =  np.cross(l_[:,ii,:],norm_)
                r_[:,:,ii,:] = polysurf[:,ii,:]-r.reshape([-1,1,3])
                R[:,:,ii] = np.sqrt(np.sum(r_[:,:,ii,:]*r_[:,:,ii,:], axis=-1))
                
            d = np.sum(r_[:,:,0,:]*norm_,axis=-1) # nr*ne
            P_ = np.zeros_like(r_) # nr*ne*3*3
            temp = d.reshape([d.shape[0],d.shape[1],1])*norm_.reshape([1,norm_.shape[0], norm_.shape[1]]) # nr*ne*3
            for ii in xrange(3):
                P_[:,:,ii,:] = r_[:,:,ii,:]-temp
            P0 = np.zeros_like(R) # nr*ne*3
            lpos = np.zeros_like(R)
            lneg = np.zeros_like(R)
            R0 = np.zeros_like(R)
            for ii in xrange(3):
                P0[:,:,ii] = np.abs(np.sum(P_[:,:,ii,:]*u_[:,ii,:],axis=-1))
                lpos[:,:,ii] = np.sum(P_[:,:,(ii+1)%3,:]*l_[:,ii,:],axis=-1)
                lneg[:,:,ii] = np.sum(P_[:,:,ii,:]*l_[:,ii,:],axis=-1)
                R0[:,:,ii] = np.sqrt(P0[:,:,ii]*P0[:,:,ii]+d*d)
            result = np.zeros([r.shape[0],polysurf.shape[0]])
            
            R0__2 = R0**2
            absd = np.abs(d)
            for ii in xrange(3):
                noise = 1.e-10*np.sqrt(np.sum(l_[:,ii,:]*l_[:,ii,:],axis=-1))                
                
                
                check2 = (R[:,:,ii]+lneg[:,:,ii])>noise
                lg = np.where(check2,\
                              np.log(R[:,:,(ii+1)%3]+lpos[:,:,ii]) - np.log(R[:,:,ii]+lneg[:,:,ii]),\
                              np.zeros([r.shape[0],polysurf.shape[0]])\
                              )
                Im1 = lg
                Ip1 = 0.5*(R[:,:,(ii+1)%3]*lpos[:,:,ii]-R[:,:,ii]*lneg[:,:,ii]+R0[:,:,ii]*R0[:,:,ii]*Im1)                
                check3 = absd>noise                
                result_branch2 = np.where(check3,\
                                          (P0[:,:,ii]*Ip1+absd**2\
                                           *(P0[:,:,ii]*Im1-absd\
                                             *(np.arctan(P0[:,:,ii]*lpos[:,:,ii]/(R0__2[:,:,ii]+absd*R[:,:,(ii+1)%3]))\
                                               -np.arctan(P0[:,:,ii]*lneg[:,:,ii]/(R0__2[:,:,ii]+absd*R[:,:,ii]))\
                                               )\
                                               )\
                                               )/3.0,\
                                          (P0[:,:,ii]*Ip1)/3.0)       
                
                check1 = (R0[:,:,ii] < noise)
                temp_result_add = np.where(check1.reshape([r.shape[0],polysurf.shape[0]]), \
                                           np.zeros_like(result), \
                                           result_branch2)   
                result = np.where( np.sum(P_[:,:,ii,:]*u_[:,ii,:], axis=-1)>0,\
                                  result + temp_result_add,\
                                  result - temp_result_add)
            return result
            pass
        except Exception as e:
            print e
            raise                


            
if __name__ == '__main__':
    x = np.arange(2)
    y = np.arange(2)
    
    import itertools
    import triangle
    import triangle.plot as triplt
    import numpy as np
    import matplotlib.pylab as plt
    grids = np.array( list(itertools.product(x,y)))
    theta = 0./180*np.pi
    rot = np.array([[np.cos(theta),np.sin(theta)],
                     [-np.sin(theta),np.cos(theta)]])
    
    grids = np.dot(grids,rot)
    
        
    r = [[0.1,0.1,0],[0.9,0.75,0.1],[1.2,0.6,0]]
    r = np.array(r)
#    print a
    plt.figure()
    ax1 = plt.subplot(111, aspect='equal')
    a = dict({'vertices':grids})
    triplt.plot(ax1,**a)
    ax1.plot(r[:,0],r[:,1],'ro')
    b = triangle.triangulate(a)
    triplt.plot(ax1,**b)
    plt.show()
    
    grid3 = np.hstack([grids, np.zeros([grids.shape[0],1])])
    trias = b['triangles']
    
    formPoly = lambda ii: grid3[trias[ii],:]
    poly = map(formPoly, xrange(trias.shape[0]))
#    print grid3[trias[3]]
#    print trias[3]
    poly = np.array(poly)
       
    sig_v = SiguralityInt_vec()
    tt1 = sig_v._1_R_Polygon(polysurf=poly, r=r)
#    tt2 = sig_v._R_Polygon(polysurf=poly, r=r)
#    tt3 = sig_v._grad_R_Polygon(polysurf=poly, r=r)
    
    sig_s = SiguralityInt()
    ftt1 = lambda cc : sig_s._1_R_Polygon(polysurf=poly[cc], r=r[0])
    ftt2 = lambda cc : sig_s._1_R_Polygon(polysurf=poly[cc], r=r[1])
    ftt3 = lambda cc : sig_s._1_R_Polygon(polysurf=poly[cc], r=r[2])
    tt1_s = map(ftt1, xrange(poly.shape[0]))
    tt2_s = map(ftt2, xrange(poly.shape[0]))
    tt3_s = map(ftt3, xrange(poly.shape[0]))
    print tt1[0,:]-np.array(tt1_s)
    print tt1[1,:]-np.array(tt2_s)
    print tt1[2,:]-np.array(tt3_s)
#    