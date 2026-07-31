<!-- Home.vue -->
<template>
  <div class="container view-container home-view">
    <div class="row py-4 px-0">
      <div class="col-12 col-lg-9 order-1 order-lg-2 typography-body mb-4 mb-lg-0" ref="mainContent">
        <div class="row">
          <div class="col">
            <div v-if="currentTag" class="mb-3 d-flex align-items-center gap-2">
              <span>{{ filteredByText }}：</span>
              <span class="current-tag-chip d-inline-flex">
                <span># {{ currentTag }}</span>
                <button class="chip-close" @click="clearTag">×</button>
              </span>
            </div>
            <PostList :docs="filteredDocs" :perPage="5" />
          </div>
        </div>
      </div>

      <div class="col-12 col-lg-3 order-2 order-lg-1" ref="sidebarContainer">
        <div class="sticky-sidebar" ref="sidebarContent">
          <div class="d-flex flex-column w-100 gap-4">
            <ProfileCard class="w-100" />
            <TagCloud class="w-100" :tagData="tagList" />
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
defineOptions({ name: 'HomeView' })
import { computed, nextTick, onBeforeUnmount, onMounted, ref, useTemplateRef, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { useHead } from '@unhead/vue'
import { useRoute, useRouter } from 'vue-router'
import ProfileCard from '@/components/layout/ProfileCard.vue'
import TagCloud from '@/components/layout/TagCloud.vue'
import PostList from '@/components/layout/PostList.vue'
import { loadPosts } from '@/utils/contentLoader'

useHead({ title: 'gblog - Home' })

const { t, locale } = useI18n()
const route = useRoute()
const router = useRouter()

const postData = ref<any[]>([])
const sidebarContent = useTemplateRef('sidebarContent')
const sidebarContainer = useTemplateRef('sidebarContainer')

const filteredByText = computed(() => t('filteredBy'))
const currentTag = computed(() => route.query.tag || '')
const filteredDocs = computed(() => {
  const tag = currentTag.value
  if (!tag) return postData.value
  return postData.value.filter(p => Array.isArray(p.tags) && p.tags.includes(tag))
})
const tagList = computed(() => {
  const set = new Set<string>()
  postData.value.forEach(p => (p.tags || []).forEach((t: any) => set.add(t)))
  return Array.from(set).sort()
})

function loadPostData() {
  try {
    postData.value = loadPosts() || []
  } catch (error) {
    console.error('Failed to load post data:', error)
    postData.value = []
  }
}

function updateSidebarDimensions() {
  if (window.innerWidth < 992) return

  const header = document.querySelector('header')
  const footer = document.querySelector('footer')
  const content = sidebarContent.value
  const sidebarContainerEl = sidebarContainer.value

  if (!content || !sidebarContainerEl) return

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

  content.style.maxHeight = `${availableHeight}px`
  content.style.overflowY = content.scrollHeight > availableHeight ? 'auto' : 'visible'
}

function clearTag() {
  const q = { ...route.query }
  delete q.tag
  q.page = '1'
  router.push({ path: route.path, query: q }).catch(() => { })
  nextTick(() => {
    window.scrollTo({ top: 0, behavior: 'smooth' })
    loadPostData()
  })
}

loadPostData()

watch(locale, (newLocale, oldLocale) => {
  if (newLocale !== oldLocale) {
    if (currentTag.value) {
      clearTag()
    } else {
      loadPostData()
    }
  }
})

onMounted(() => {
  updateSidebarDimensions()
  window.addEventListener('scroll', updateSidebarDimensions)
  window.addEventListener('resize', updateSidebarDimensions)
})

onBeforeUnmount(() => {
  window.removeEventListener('scroll', updateSidebarDimensions)
  window.removeEventListener('resize', updateSidebarDimensions)
})
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

.current-tag-chip {
  font-size: 1rem;
  font-weight: 500;
  color: var(--app-chip-text);
  background: var(--app-chip-bg);
  padding: 0.25rem 0.5rem;
  border: 1px solid var(--app-chip-border);
  box-shadow: var(--app-card-shadow);
  border-radius: 12px;
  align-items: center;
  gap: 0.4rem;
}

.chip-close {
  font-size: 1.2rem;
  line-height: 1;
  background: transparent;
  border: none;
  color: var(--app-chip-close-text);
  padding: 0;
  margin-left: 2px;
  cursor: pointer;
  opacity: 1;
  transition: color 0.2s ease;
}

.chip-close:hover {
  color: var(--app-primary);
}
</style>
