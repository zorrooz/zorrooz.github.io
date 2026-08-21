import HomeView from '@/views/Home.vue'
import CategoryView from '@/views/Category.vue'
import ResourceView from '@/views/Resource.vue'
import AboutView from '@/views/About.vue'
import ArticleView from '@/views/Article.vue'
import NotFoundView from '@/views/NotFound.vue'

import type { RouteRecordRaw } from 'vue-router'

import { LOCALE_SEGMENTS, preferredLocaleSegment, type LocaleSegment } from '@/config'

const prefixedRoutes = (prefix: LocaleSegment): RouteRecordRaw[] => [
  { path: `/${prefix}/`, name: `${prefix}-Home`, component: HomeView },
  { path: `/${prefix}/category`, name: `${prefix}-Category`, component: CategoryView },
  { path: `/${prefix}/resource`, name: `${prefix}-Resource`, component: ResourceView },
  { path: `/${prefix}/about`, name: `${prefix}-About`, component: AboutView },
  {
    path: `/${prefix}/article/:path*`,
    name: `${prefix}-Article`,
    component: ArticleView,
    props: true,
  },
]

export const routes: RouteRecordRaw[] = [
  // 无前缀旧 URL 重定向到带 locale 前缀的等价路径
  { path: '/', redirect: () => `/${preferredLocaleSegment()}/` },
  { path: '/zh', redirect: '/zh/' },
  { path: '/en', redirect: '/en/' },
  { path: '/category', redirect: (to) => `/${preferredLocaleSegment()}${to.path}` },
  { path: '/resource', redirect: (to) => `/${preferredLocaleSegment()}${to.path}` },
  { path: '/about', redirect: (to) => `/${preferredLocaleSegment()}${to.path}` },
  { path: '/article/:path*', redirect: (to) => `/${preferredLocaleSegment()}${to.path}` },
  // 其余无前缀路径：先补 locale 前缀，落入下方 404
  { path: '/:pathMatch(.*)*', redirect: (to) => `/${preferredLocaleSegment()}${to.path}` },
  ...LOCALE_SEGMENTS.map(prefixedRoutes).flat(),
  // 语言前缀下的未知路径 → 品牌 404
  { path: '/:locale(zh|en)/:pathMatch(.*)*', name: 'NotFound', component: NotFoundView },
]
