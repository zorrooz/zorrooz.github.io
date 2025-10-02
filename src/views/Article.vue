<template>
  <div class="container">
    <div class="col-12 col-lg-2">
      <!-- <NavigationTree /> -->
    </div>
    <div class="row justify-content-center px-3">
      <div class="col-12 col-lg-8">
        <RenderMarkdown :rawMarkdown="rawMarkdown" />
      </div>

      <div class="col-12 col-lg-2" ref="sidebarContainer">
        <!-- 使用通用的顶部间距类 -->
        <div class="sticky-sidebar sidebar-top-spacing" ref="stickySidebar">
          <div class="sidebar-content" ref="sidebarContent">
            <MarkdownToc containerSelector=".markdown-body" :levels="[2, 3, 4, 5, 6]" />
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import RenderMarkdown from '@/components/layout/RenderMarkdown.vue'
import rawMarkdown from '@/content/demo.md?raw'
import MarkdownToc from '@/components/layout/MarkdownToc.vue'
// import NavigationTree from '@/components/layout/NavigationTree.vue'

export default {
  name: 'ArticleView',
  components: { RenderMarkdown, MarkdownToc },
  data() {
    return {
      rawMarkdown,
      // 通用命名：侧边栏顶部间距（像素），可根据需求调整
      sidebarTopSpacing: 8,
      isScrolled: false
    }
  },
  mounted() {
    // 初始化粘性定位的顶部距离
    this.setStickyPosition()
    this.updateSidebarLayout()

    window.addEventListener('scroll', this.handleScroll)
    window.addEventListener('resize', this.updateSidebarLayout)
  },
  beforeUnmount() {
    window.removeEventListener('scroll', this.handleScroll)
    window.removeEventListener('resize', this.updateSidebarLayout)
  },
  methods: {
    // 设置粘性定位参数
    setStickyPosition() {
      const stickySidebar = this.$refs.stickySidebar
      if (stickySidebar) {
        // 粘性定位时的顶部距离与初始间距保持一致
        stickySidebar.style.top = `${this.sidebarTopSpacing}px`
      }
    },

    // 处理滚动事件
    handleScroll() {
      this.isScrolled = window.scrollY > 0
      this.updateSidebarLayout()
    },

    // 更新侧边栏布局（更通用的方法名）
    updateSidebarLayout() {
      const header = document.querySelector('header')
      const footer = document.querySelector('footer')
      const sidebarContent = this.$refs.sidebarContent
      const sidebarContainer = this.$refs.sidebarContainer
      const stickySidebar = this.$refs.stickySidebar

      if (!sidebarContent || !sidebarContainer || !stickySidebar) return

      const headerHeight = header?.offsetHeight || 0
      const footerHeight = footer?.offsetHeight || 0
      const viewportHeight = window.innerHeight
      const scrollTop = window.scrollY
      const documentHeight = document.documentElement.scrollHeight

      // 计算剩余页面高度（考虑顶部间距）
      const remainingPageHeight = Math.max(
        0,
        documentHeight - scrollTop - headerHeight - footerHeight - 40 - this.sidebarTopSpacing
      )

      // 计算可用高度（考虑顶部间距）
      const availableHeight = Math.min(
        viewportHeight - headerHeight - 40 - this.sidebarTopSpacing,
        remainingPageHeight
      )

      // 设置内容区高度和滚动状态
      sidebarContent.style.maxHeight = `${availableHeight}px`
      sidebarContent.style.overflowY = sidebarContent.scrollHeight > availableHeight ? 'auto' : 'visible'

      // 计算底部定位
      const containerTop = sidebarContainer.getBoundingClientRect().top + scrollTop
      const maxStickyBottom = documentHeight - footerHeight - 20
      const sidebarBottom = containerTop + availableHeight + this.sidebarTopSpacing

      stickySidebar.style.bottom = sidebarBottom > maxStickyBottom
        ? `${sidebarBottom - maxStickyBottom}px`
        : '20px'
    }
  }
}
</script>

<style scoped>
.sticky-sidebar {
  position: sticky;
  box-sizing: border-box;
  width: 100%;
  transition: bottom 0.2s ease;
}

/* 通用的顶部间距样式，不与特定类名绑定 */
.sidebar-top-spacing {
  margin-top: 0.5rem !important;
  /* 可根据实际需求调整 */
}

.sidebar-content {
  width: 100%;
  box-sizing: border-box;
  padding: 15px;
  background: #fff;
  border-radius: 8px;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.05);
  -webkit-overflow-scrolling: touch;
  transition: max-height 0.2s ease;
}

@media (max-width: 991px) {
  .sticky-sidebar {
    position: static;
    top: auto !important;
    bottom: auto !important;
  }

  .sidebar-content {
    max-height: none !important;
    overflow-y: visible !important;
  }
}
</style>
