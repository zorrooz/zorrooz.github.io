<!-- Home.vue -->
<template>
  <div class="container">
    <div class="row pt-4 px-0">
      <!-- 中间主内容区域 -->
      <div class="col-12 col-lg-9 order-1 order-lg-2 typography-body" ref="mainContent">
        <div class="row">
          <div class="col">
            <p class="mb-2">abcd 测试</p>
            <h3 class="mb-4">文档中心</h3>
            <PostList :docs="postData" :perPage="5" />
          </div>
        </div>
      </div>

      <!-- 左侧边栏 -->
      <div class="col-12 col-lg-3 order-2 order-lg-1 pb-2 px-2" ref="sidebarContainer">
        <div class="sticky-sidebar" ref="sidebarContent">
          <div class="d-flex flex-column align-items-center align-items-lg-end w-100 gap-3">
            <ProfileCard class="w-100" />
            <TagCloud class="w-100" />
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import ProfileCard from '@/components/layout/ProfileCard.vue'
import TagCloud from '@/components/layout/TagCloud.vue'
import PostList from '@/components/layout/PostList.vue'
import postData from '@/content/posts.json'

export default {
  name: 'HomeView',
  components: { ProfileCard, TagCloud, PostList },
  data() { return { postData } },
  mounted() {
    this.updateSidebarDimensions()
    window.addEventListener('scroll', this.updateSidebarDimensions)
    window.addEventListener('resize', this.updateSidebarDimensions)
  },
  beforeUnmount() {
    window.removeEventListener('scroll', this.updateSidebarDimensions)
    window.removeEventListener('resize', this.updateSidebarDimensions)
  },
  methods: {
    updateSidebarDimensions() {
      const header = document.querySelector('header')
      const footer = document.querySelector('footer')
      const content = this.$refs.sidebarContent
      const sidebarContainer = this.$refs.sidebarContainer

      if (!content || !sidebarContainer) return

      const headerHeight = header?.offsetHeight || 0
      const footerHeight = footer?.offsetHeight || 0
      const viewportHeight = window.innerHeight
      const scrollTop = window.scrollY
      const documentHeight = document.documentElement.scrollHeight

      const remainingPageHeight = Math.max(
        0,
        documentHeight - scrollTop - headerHeight - footerHeight - 40
      )
      const availableHeight = Math.min(
        viewportHeight - headerHeight - 40,
        remainingPageHeight
      )

      const containerTop = sidebarContainer.getBoundingClientRect().top + scrollTop
      const maxStickyBottom = documentHeight - footerHeight - 20

      content.style.maxHeight = `${availableHeight}px`
      content.style.overflowY = content.scrollHeight > availableHeight ? 'auto' : 'visible'

      const sidebar = content.parentElement
      if (sidebar) {
        const sidebarBottom = containerTop + availableHeight
        sidebar.style.bottom = sidebarBottom > maxStickyBottom
          ? `${sidebarBottom - maxStickyBottom}px`
          : ''
      }
    }
  }
}
</script>

<style scoped>
.sticky-sidebar {
  position: sticky;
  top: 30px;
  box-sizing: border-box;
  width: 100%;
  -webkit-overflow-scrolling: touch;
  transition: max-height 0.2s ease;
}

@media (max-width: 991px) {
  .sticky-sidebar {
    position: static;
    top: auto;
    bottom: auto !important;
    max-height: none !important;
    overflow-y: visible !important;
  }
}
</style>
