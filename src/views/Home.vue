<!-- Home.vue -->
<template>
  <div class="container">
    <div class="row py-4 px-0">
      <!-- 中间主内容区域 -->
      <div class="col-12 col-lg-9 order-1 order-lg-2 typography-body" ref="mainContent">
        <div class="row">
          <div class="col">


            <div v-if="currentTag" class="mb-3">
              <span class="current-tag-chip d-inline-flex">
                <span>标签：# {{ currentTag }}</span>
                <button class="chip-close" @click="clearTag" aria-label="清除筛选" title="清除筛选">×</button>
              </span>
            </div>
            <PostList :docs="filteredDocs" :perPage="10" />
          </div>
        </div>
      </div>

      <!-- 左侧边栏 -->
      <div class="col-12 col-lg-3 order-2 order-lg-1 pb-2" ref="sidebarContainer">
        <div class="row">
          <div class="col">
            <div class="sticky-sidebar" ref="sidebarContent">
              <div class="d-flex flex-column align-items-center align-items-lg-end w-100 gap-0">
                <ProfileCard class="w-100" />
                <TagCloud class="w-100" :tagData="tagList" />
              </div>
            </div>
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
  computed: {
    currentTag() { return this.$route.query.tag || '' },
    filteredDocs() {
      const tag = this.currentTag
      if (!tag) return this.postData
      return this.postData.filter(p => Array.isArray(p.tags) && p.tags.includes(tag))
    },
    tagList() {
      const set = new Set()
      this.postData.forEach(p => (p.tags || []).forEach(t => set.add(t)))
      return Array.from(set).sort()
    }
  },
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
    },
    clearTag() {
      const q = { ...this.$route.query }
      delete q.tag
      q.page = '1'
      this.$router.push({ path: this.$route.path, query: q }).catch(() => {})
      this.$nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }))
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
.current-tag-chip {
  font-size: 1.35rem;
  font-weight: 400;
  color: #000; /* 纯黑色显示 */
  background: transparent;
  padding: 0;
  border: none;
  box-shadow: none;
  border-radius: 0;
  align-items: baseline; /* 与文字基线对齐 */
  gap: 0.2rem; /* 更紧凑的间距 */
}

.chip-close {
  font-size: 1.6rem; /* 更大一点，与文字同水平更醒目 */
  line-height: 1;
  background: transparent;
  border: none;
  color: #6c757d; /* 与全局 text-secondary 接近 */
  padding: 0;
  cursor: pointer;
  opacity: 0;
  transition: opacity 0.2s ease;
}

.current-tag-chip:hover .chip-close {
  opacity: 1;
}

.chip-close:hover {
  color: #047AFF; /* 与全局主色一致的悬停色 */
}
</style>
