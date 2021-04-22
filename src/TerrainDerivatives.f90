!-------------------------------------------------------------------

subroutine slope (h,s,nx,ny,dx,dy)

  implicit none

  integer nx,ny
  double precision h(nx*ny),s(nx*ny),dx,dy

  integer i,j,ij,ia,ib,ic,id,ie,if,ig,ih,ii
  double precision dzdx,dzdy,con

  con=45.d0/atan(1.d0)

  s=0.d0
  do j=2,ny-1
    do i=2,nx-1
      ij=i+(j-1)*nx
      ia=ij+nx-1
      ib=ia+1
      ic=ib+1
      id=ij-1
      ie=ij
      if=ij+1
      ig=ij-nx-1
      ih=ig+1
      ii=ih+1
      dzdx=((h(ic)+2.d0*h(if)+h(ii))-(h(ia)+2.d0*h(id)+h(ig)))/8.d0/dx
      dzdy=((h(ig)+2.d0*h(ih)+h(ii))-(h(ia)+2.d0*h(ib)+h(ic)))/8.d0/dy
      s(ij)=dzdx**2+dzdy**2
      if (s(ij).gt.tiny(s(ij))) s(ij)=atan(sqrt(s(ij)))*con
    enddo
  enddo

end subroutine slope

!-------------------------------------------------------------------

subroutine curvature (h,curv,nx,ny,dx,dy)

  implicit none

  integer nx,ny
  double precision dx,dy
  double precision h(nx*ny),curv(nx*ny)

  integer i,j,ij,i1,i2,i3,i4,i5,i6,i7,i8,i9
  double precision a,b,c,d,e,f

  curv=0.d0
  do j=2,ny-1
    do i=2,nx-1
      ij=i+(j-1)*nx
      i1=ij+nx-1
      i2=i1+1
      i3=i2+1
      i4=ij-1
      i5=ij
      i6=ij+1
      i7=ij-nx-1
      i8=i7+1
      i9=i8+1
      a=(h(i1)+h(i3)+h(i4)+h(i6)+h(i7)+h(i9))/dx/dx/12.d0-(h(i2)+h(i5)+h(i8))/dx/dx/6.d0
      b=(h(i1)+h(i2)+h(i3)+h(i7)+h(i8)+h(i9))/dy/dy/12.d0-(h(i4)+h(i5)+h(i6))/dy/dy/6.d0
      c=(h(i3)+h(i7)-h(i1)-h(i9))/dx/dy/4.d0
      d=(h(i3)+h(i6)+h(i9)-h(i1)-h(i4)-h(i7))/dx/6.d0
      e=(h(i1)+h(i2)+h(i3)-h(i7)-h(i8)-h(i9))/dy/6.d0
      f=(2.d0*(h(i2)+h(i4)+h(i6)+h(i8))-(h(i1)+h(i3)+h(i7)+h(i9))+5.d0*h(i5))/9.d0
      curv(ij)=1.d0+d**2+e**2
      if (curv(ij).gt.tiny(curv(ij))) curv(ij)=(a*(1.d0+e**2)+b*(1.d0+d**2)-c*d*e)/(curv(ij)**(3.d0/2.d0))
    enddo
  enddo

end subroutine curvature
